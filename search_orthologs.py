#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_orthologs.py — Ortholog gene search and classification using miniprot.

Aligns OrthoDB protein profiles against a genome assembly with miniprot
(Heng Li), then classifies each ortholog group as Complete, Fragmented,
or Missing based on alignment coverage and score thresholds.

Supports OrthoDB lineage datasets from ODB9, ODB10, and ODB11.

Output: per-contig TSV files consumable by hapsolo.py's importBuscos().

Dependencies: miniprot (https://github.com/lh3/miniprot)
"""
import argparse
import glob
import os
import subprocess
import sys


def find_protein_file(lineage_dir):
    """Locate the protein sequence FASTA in an OrthoDB lineage dataset.

    Tries ODB10/11 format first, then ODB9 format.
    Returns the path to the protein FASTA file.
    """
    candidates = [
        os.path.join(lineage_dir, 'refseq_db.faa'),          # ODB10/11
        os.path.join(lineage_dir, 'ancestral_variants'),      # ODB9 (multi-seq)
        os.path.join(lineage_dir, 'ancestral'),               # ODB9 (consensus)
    ]
    for path in candidates:
        if os.path.exists(path):
            return path

    # Try glob for any .faa file
    faa_files = glob.glob(os.path.join(lineage_dir, '*.faa'))
    if faa_files:
        return faa_files[0]

    return None


def load_scores_cutoff(lineage_dir):
    """Load the per-BUSCO score thresholds from the lineage dataset.

    Returns dict: {busco_id: min_score}
    """
    path = os.path.join(lineage_dir, 'scores_cutoff')
    if not os.path.exists(path):
        return {}

    cutoffs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) >= 2:
                busco_id = fields[0]
                try:
                    score = float(fields[1])
                except ValueError:
                    continue
                cutoffs[busco_id] = score
    return cutoffs


def load_lengths_cutoff(lineage_dir):
    """Load the per-BUSCO expected protein lengths from the lineage dataset.

    Handles both ODB9 and ODB10/11 column orders.
    Returns dict: {busco_id: (mean_length, std_dev)}
    """
    path = os.path.join(lineage_dir, 'lengths_cutoff')
    if not os.path.exists(path):
        return {}

    cutoffs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) >= 3:
                busco_id = fields[0]
                # Try to determine column order by checking which fields are numeric
                # ODB9: BUSCO_ID  mean  sd  n_species
                # ODB10: BUSCO_ID  n_species  mean  sd
                try:
                    val1 = float(fields[1])
                    val2 = float(fields[2])
                    if len(fields) >= 4:
                        val3 = float(fields[3])
                        # If field[1] is small integer (n_species) and field[2] is larger (mean length)
                        if val1 == int(val1) and val1 < 1000 and val2 > val1:
                            # ODB10 format: ID, n_species, mean, sd
                            cutoffs[busco_id] = (val2, val3)
                        else:
                            # ODB9 format: ID, mean, sd, n_species
                            cutoffs[busco_id] = (val1, val2)
                    else:
                        # 3 columns: ID, mean, sd
                        cutoffs[busco_id] = (val1, val2)
                except ValueError:
                    continue
    return cutoffs


def get_busco_ids_from_proteins(protein_file):
    """Extract BUSCO IDs from the protein FASTA headers.

    Returns set of BUSCO IDs (the part before the colon or first underscore
    in the FASTA header, depending on format).
    """
    ids = set()
    with open(protein_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                # ODB10: headers like ">100at7147:protein_name"
                # ODB9: headers like ">EOG09150001"
                busco_id = header.split(':')[0]
                ids.add(busco_id)
    return ids


def run_miniprot(genome_fasta, protein_fasta, output_paf, threads=1):
    """Run miniprot to align proteins against the genome.

    Returns the path to the PAF output file.
    """
    cmd = [
        'miniprot',
        '-t', str(threads),
        '--paf',
        '-I',   # no secondary alignments (report best only)
        genome_fasta,
        protein_fasta,
    ]

    print('Running: ' + ' '.join(cmd))
    with open(output_paf, 'w') as fout:
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE)

    if proc.returncode != 0:
        print('miniprot failed:')
        print(proc.stderr.decode('utf-8', errors='replace'))
        sys.exit(1)

    return output_paf


def parse_miniprot_paf(paf_file):
    """Parse miniprot PAF output into per-BUSCO, per-contig hits.

    Returns list of dicts: [{busco_id, contig, query_len, query_start, query_end,
                             target_start, target_end, score, aligned_len}, ...]
    """
    hits = []
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            query_name = fields[0]
            # Extract BUSCO ID from protein name (before colon or full name)
            busco_id = query_name.split(':')[0]

            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            contig = fields[5]
            target_start = int(fields[7])
            target_end = int(fields[8])

            # Get alignment score from AS tag
            score = 0
            for tag in fields[12:]:
                if tag.startswith('AS:i:'):
                    score = int(tag[5:])
                    break

            aligned_len = query_end - query_start  # in protein coordinates

            hits.append({
                'busco_id': busco_id,
                'contig': contig,
                'query_len': query_len,
                'query_start': query_start,
                'query_end': query_end,
                'target_start': target_start,
                'target_end': target_end,
                'score': score,
                'aligned_len': aligned_len,
            })

    return hits


def classify_buscos(hits, all_busco_ids, scores_cutoff, lengths_cutoff):
    """Classify BUSCO hits as Complete, Fragmented, or Missing per contig.

    Returns dict: {contig: {busco_id: (status, start, end, score, length)}}
    """
    # Group hits by BUSCO ID
    hits_by_busco = {}
    for hit in hits:
        bid = hit['busco_id']
        if bid not in hits_by_busco:
            hits_by_busco[bid] = []
        hits_by_busco[bid].append(hit)

    # Classify each BUSCO
    # result[contig][busco_id] = (status, start, end, score, length)
    contig_results = {}

    for busco_id in all_busco_ids:
        if busco_id not in hits_by_busco:
            # No hits at all -> Missing (will be added per-contig later)
            continue

        busco_hits = hits_by_busco[busco_id]

        # Get classification thresholds
        min_score = scores_cutoff.get(busco_id, 0)
        if busco_id in lengths_cutoff:
            mean_len, std_dev = lengths_cutoff[busco_id]
            complete_threshold = mean_len - 2 * std_dev
        else:
            # No length info: use 95% of query length from best hit
            best_qlen = max(h['query_len'] for h in busco_hits)
            complete_threshold = best_qlen * 0.95

        # Filter by score cutoff
        significant_hits = [h for h in busco_hits if h['score'] >= min_score]

        if not significant_hits:
            # No significant hits -> Missing
            continue

        for hit in significant_hits:
            contig = hit['contig']
            if contig not in contig_results:
                contig_results[contig] = {}

            aligned_len = hit['aligned_len']

            if aligned_len >= complete_threshold:
                status = 'Complete'
            else:
                status = 'Fragmented'

            # Keep the best hit per BUSCO per contig
            if busco_id not in contig_results[contig]:
                contig_results[contig][busco_id] = (
                    status, hit['target_start'], hit['target_end'],
                    hit['score'], aligned_len)
            else:
                existing = contig_results[contig][busco_id]
                if hit['score'] > existing[3]:
                    contig_results[contig][busco_id] = (
                        status, hit['target_start'], hit['target_end'],
                        hit['score'], aligned_len)
                # Upgrade Fragmented to Complete if better hit found
                elif status == 'Complete' and existing[0] == 'Fragmented':
                    contig_results[contig][busco_id] = (
                        status, hit['target_start'], hit['target_end'],
                        hit['score'], aligned_len)

    return contig_results


def write_busco_output(contig_results, all_busco_ids, output_dir, lineage_name,
                       contig_fasta_dir='contigs'):
    """Write per-contig BUSCO V3 format TSV files.

    Creates the directory structure expected by hapsolo.py's importBuscos():
      output_dir/busco_CONTIG/run_CONTIG/full_table_CONTIG.tsv
    """
    all_contigs = set(contig_results.keys())

    for contig in sorted(all_contigs):
        busco_dir = os.path.join(output_dir, 'busco_' + contig, 'run_' + contig)
        os.makedirs(busco_dir, exist_ok=True)

        tsv_path = os.path.join(busco_dir, 'full_table_' + contig + '.tsv')
        with open(tsv_path, 'w') as f:
            f.write('# BUSCO version is: search_orthologs 1.0 (miniprot)\n')
            f.write('# The lineage dataset is: ' + lineage_name + '\n')
            f.write('# To reproduce this run: python search_orthologs.py -i '
                    + contig_fasta_dir + '/' + contig + '.fasta -l '
                    + lineage_name + '\n')
            f.write('#\n')
            f.write('# Busco id\tStatus\tContig\tStart\tEnd\tScore\tLength\n')

            contig_buscos = contig_results.get(contig, {})
            for busco_id in sorted(all_busco_ids):
                if busco_id in contig_buscos:
                    status, start, end, score, length = contig_buscos[busco_id]
                    f.write(busco_id + '\t' + status + '\t' + contig + '\t'
                            + str(start) + '\t' + str(end) + '\t'
                            + str(score) + '\t' + str(length) + '\n')
                else:
                    f.write(busco_id + '\tMissing\n')

    print('BUSCO results written for ' + str(len(all_contigs)) + ' contigs to ' + output_dir)


def detect_lineage_name(lineage_dir):
    """Extract the lineage name from the directory path."""
    name = os.path.basename(os.path.normpath(lineage_dir))
    # Check for dataset.cfg
    cfg = os.path.join(lineage_dir, 'dataset.cfg')
    if os.path.exists(cfg):
        with open(cfg) as f:
            for line in f:
                if line.startswith('name'):
                    parts = line.strip().split('=')
                    if len(parts) >= 2:
                        return parts[1].strip()
    return name


def main():
    parser = argparse.ArgumentParser(
        description='Lightweight BUSCO classifier using miniprot. '
                    'Replaces BUSCO for HapSolo pipeline. '
                    'Supports ODB9, ODB10, and ODB11 lineage datasets.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input genome FASTA file')
    parser.add_argument('-l', '--lineage', required=True,
                        help='Path to OrthoDB lineage dataset directory '
                             '(e.g., diptera_odb10/)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for BUSCO results')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads for miniprot (default: 1)')
    parser.add_argument('--contig-dir', default='contigs',
                        help='Contig FASTA directory name for output headers '
                             '(default: contigs)')
    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.input):
        print('Error: Input file not found: ' + args.input)
        sys.exit(1)
    if not os.path.isdir(args.lineage):
        print('Error: Lineage directory not found: ' + args.lineage)
        sys.exit(1)

    # Find protein sequences
    protein_file = find_protein_file(args.lineage)
    if protein_file is None:
        print('Error: No protein sequence file found in lineage directory.')
        print('Expected one of: refseq_db.faa, ancestral_variants, ancestral')
        sys.exit(1)
    print('Protein sequences: ' + protein_file)

    # Load classification thresholds
    scores_cutoff = load_scores_cutoff(args.lineage)
    lengths_cutoff = load_lengths_cutoff(args.lineage)
    print('Score cutoffs loaded: ' + str(len(scores_cutoff)) + ' BUSCOs')
    print('Length cutoffs loaded: ' + str(len(lengths_cutoff)) + ' BUSCOs')

    # Get all BUSCO IDs
    all_busco_ids = get_busco_ids_from_proteins(protein_file)
    print('Total BUSCO genes in lineage: ' + str(len(all_busco_ids)))

    # Detect lineage name
    lineage_name = detect_lineage_name(args.lineage)
    print('Lineage: ' + lineage_name)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Run miniprot alignment
    paf_file = os.path.join(args.output, 'miniprot_busco.paf')
    run_miniprot(args.input, protein_file, paf_file, args.threads)

    # Parse results
    hits = parse_miniprot_paf(paf_file)
    print('Total alignment hits: ' + str(len(hits)))

    # Classify
    contig_results = classify_buscos(hits, all_busco_ids, scores_cutoff,
                                     lengths_cutoff)

    # Write output in BUSCO V3 format
    write_busco_output(contig_results, all_busco_ids, args.output,
                       lineage_name, args.contig_dir)

    # Print summary
    total_complete = 0
    total_fragmented = 0
    total_missing = 0
    all_contigs_buscos = {}
    for contig, buscos in contig_results.items():
        for bid, (status, s, e, sc, l) in buscos.items():
            if bid not in all_contigs_buscos or status == 'Complete':
                all_contigs_buscos[bid] = status
    for bid in all_busco_ids:
        st = all_contigs_buscos.get(bid, 'Missing')
        if st == 'Complete':
            total_complete += 1
        elif st == 'Fragmented':
            total_fragmented += 1
        else:
            total_missing += 1

    print('\nSummary:')
    print('  Complete:    ' + str(total_complete))
    print('  Fragmented:  ' + str(total_fragmented))
    print('  Missing:     ' + str(total_missing))
    print('  Total:       ' + str(len(all_busco_ids)))


if __name__ == '__main__':
    main()
