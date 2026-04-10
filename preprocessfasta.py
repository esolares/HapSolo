#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import re, os, argparse

# usage preprocessfasta.py -i myfastafile.fasta
parser = argparse.ArgumentParser(description='Preprocess FASTA file and outputs a clean FASTA and seperates contigs based on unique headers. Removes special chars')
parser.add_argument('-i', '--input', help='Input FASTA file', type=str, required=True)
parser.add_argument('-m', '--maxcontig', help='Max size of contig in Mb to output in contigs folder.', type=str, required=False)

args = parser.parse_args()
asmfilename = args.input
maxcontigsize = args.maxcontig

if maxcontigsize is None:
    maxcontigsize = 10*1000000
else:
    maxcontigsize = int(maxcontigsize)*1000000

# First pass: collect original headers (streaming, no full file load)
original_headers = []
with open(asmfilename) as fin:
    for line in fin:
        line = line.strip()
        if line and line[0] == '>':
            original_headers.append(line[1:].split()[0])

if len(original_headers) == 0:
    print('No sequences found in FASTA file: ' + asmfilename)
    quit(1)

# Sanitize headers: replace special characters with underscores
sanitized_headers = [re.sub('[^a-zA-Z0-9.]', '_', h) for h in original_headers]

# Find minimum prefix length that keeps all headers unique
maxlen = max(len(h) for h in sanitized_headers)
uniqueheadersize = -1
for i in range(1, maxlen+1):
    prefixes = set(h[0:i] for h in sanitized_headers)
    if len(prefixes) == len(sanitized_headers):
        uniqueheadersize = i
        break

if uniqueheadersize == -1:
    print('This FASTA file does not contain unique headers. Please fix and rerun again.')
    quit(1)

# Truncate to unique prefix
final_names = [h[0:uniqueheadersize] for h in sanitized_headers]

# If first underscore-delimited field alone is unique, use that instead
first_fields = set(h.split('_')[0] for h in final_names)
if len(first_fields) == len(final_names):
    final_names = [h.split('_')[0] for h in final_names]

# Build lookup: original header -> final sanitized name
name_lookup = dict()
for i in range(len(original_headers)):
    name_lookup[original_headers[i]] = final_names[i]

# Prepare output
fileext = asmfilename.split('.')[-1]
outfile = asmfilename.replace('.' + fileext, '') + '_new.' + fileext
outdir = "contigs"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Second pass: stream through FASTA and write outputs (one contig at a time)
contig_idx = 0
current_name = None
seq_parts = []

def write_contig(name, seq_parts, fout, outdir, maxcontigsize):
    seq = ''.join(seq_parts)
    fout.write('>' + name + '\n')
    fout.write(seq + '\n')
    if len(seq) < maxcontigsize:
        with open(os.path.join(outdir, name + '.fasta'), 'w') as contigout:
            contigout.write('>' + name + '\n')
            contigout.write(seq + '\n')

with open(outfile, 'w') as fout:
    with open(asmfilename) as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                # Write previous contig if any
                if current_name is not None:
                    write_contig(current_name, seq_parts, fout, outdir, maxcontigsize)
                # Start new contig
                original = line[1:].split()[0]
                current_name = name_lookup[original]
                seq_parts = []
            else:
                seq_parts.append(line)

    # Write last contig
    if current_name is not None:
        write_contig(current_name, seq_parts, fout, outdir, maxcontigsize)

# Write name mapping file: original_name -> sanitized_name
with open(os.path.join(outdir, 'name_mapping.tsv'), 'w') as mapout:
    mapout.write('#original_name\tsanitized_name\n')
    for i in range(len(original_headers)):
        mapout.write(original_headers[i] + '\t' + final_names[i] + '\n')
print('Name mapping written to ' + outdir + '/name_mapping.tsv')
