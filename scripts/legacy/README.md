# Legacy HapSolo Scripts

These scripts use the **original** HapSolo workflow:
- BLAT all-by-all alignment via array jobs (one job per contig)
- BUSCO V3 with Augustus for ortholog classification
- ODB9 lineage datasets

They are kept here for reproducing published results. New users should use the modern pipeline scripts in the parent directory (`scripts/sbatch_*.sh`, `scripts/qsub_*.sh`), which use minimap2 and miniprot directly with no BUSCO/Augustus dependency.

## Files

| Script | Purpose |
|---|---|
| `sbatch_blat.sh` / `qsub_blat.sh` | BLAT array job (one task per contig) |
| `bash_andreaconcatpsl.sh` | Concatenate per-task PSL files into a single PSL |
| `bash_gnuparallelblat.sh` | GNU Parallel version of the BLAT step |
| `sbatch_busco.sh` | BUSCO V3 array job (one task per contig) |
| `bash_buscopreprocess.sh` | Build the BUSCO job file list |
| `bash_gnuparallelbusco.sh` | GNU Parallel version of the BUSCO step |
| `bash_quastbusco.sh` | Per-contig QUAST + BUSCO worker invoked by the parallel script |

## Workflow (legacy)

1. Run `preprocessfasta.py` (or `bash_buscopreprocess.sh`) to split contigs.
2. Run `bash_buscopreprocess.sh` to create `contigs/buscojobfile.txt`.
3. Submit `sbatch_blat.sh` (or `qsub_blat.sh`) as an array job, one task per contig.
4. After BLAT completes, run `bash_andreaconcatpsl.sh` to concatenate the PSL files.
5. Submit `sbatch_busco.sh` as an array job for the BUSCO V3 classification.
6. Run `hapsolo.py` with the resulting PSL and the `contigs/busco/` directory.

## Why this is no longer recommended

- **BLAT array jobs** scale poorly: one job per contig means thousands of small jobs, lots of scheduler overhead, and inconsistent runtimes. minimap2 in a single multi-threaded job is faster end-to-end.
- **BUSCO V3** requires Python 2 and Augustus, which is hard to install on modern systems. The Singularity image (`hapsolo_busco3:0.01`) is the only practical way to use it.
- **search_orthologs.py** with miniprot replaces BUSCO V3 entirely while producing output in the same TSV format that `hapsolo.py` consumes.
