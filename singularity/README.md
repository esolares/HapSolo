# HapSolo Singularity Image

This directory contains the recipe and binary placeholders for building a self-contained HapSolo Singularity image.

## Directory Layout

```
singularity/
├── hapsolo.def          Singularity recipe
├── build.sh             Helper build script
├── README.md            This file
└── binaries/
    ├── minimap2/        Required — place minimap2 binary here
    ├── miniprot/        Required — place miniprot binary here
    └── blat/            Optional — alternative aligner, place blat binary here
```

Each subdirectory under `binaries/` has its own `README.md` with instructions for downloading or compiling that tool.

## Build Steps

1. **Place the required binaries.** At minimum you need `minimap2` and `miniprot`:

   ```
   # See singularity/binaries/minimap2/README.md
   cp /path/to/minimap2 singularity/binaries/minimap2/

   # See singularity/binaries/miniprot/README.md
   cp /path/to/miniprot singularity/binaries/miniprot/
   ```

   The optional `blat` aligner needs to be uncommented in `hapsolo.def` if you want it in the image.

2. **Build the image.** Singularity 3.x requires root for `build`:

   ```
   cd singularity/
   ./build.sh                       # uses default name hapsolo.sif
   ./build.sh hapsolo_v2.0.sif      # custom output name
   ```

   This produces `hapsolo.sif` (or your chosen name) — a single self-contained file you can copy to any HPC node with Singularity installed.

## Inside the Image

| Path | Contents |
|---|---|
| `/opt/HapSolo/` | All Python scripts (`hapsolo.py`, `hapsolo_cli.py`, `preprocessfasta.py`, `search_orthologs.py`) |
| `/opt/minimap2/minimap2` | minimap2 binary |
| `/opt/miniprot/miniprot` | miniprot binary |
| `/opt/blat/blat` | (optional) blat binary |
| `/usr/local/bin/hapsolo` | Symlink to `hapsolo_cli.py` for easy invocation |

All binaries are added to `$PATH` via the `%environment` block.

## Usage

The image's default runscript is `hapsolo_cli.py`:

```
# Bind your data directory and run the CLI directly
singularity run --bind /data:/data hapsolo.sif preprocess -i /data/assembly.fasta
singularity run --bind /data:/data hapsolo.sif align      -i /data/assembly_new.fasta -t 8
singularity run --bind /data:/data hapsolo.sif search     -i /data/assembly_new.fasta -l /data/diptera_odb10/ -t 8
singularity run --bind /data:/data hapsolo.sif train      -i /data/assembly_new.fasta --paf /data/assembly_new_self_align.paf.gz -b /data/ortholog_output/ -t 32 -n 1000
```

To run any individual tool:

```
singularity exec hapsolo.sif minimap2 --version
singularity exec hapsolo.sif miniprot --version
singularity exec hapsolo.sif hapsolo.py --help
```

To open a shell inside the image (useful for debugging):

```
singularity shell --bind /data:/data hapsolo.sif
```

## Notes

- Use `--bind /host/path:/container/path` to make host directories visible inside the container. Without `--bind`, only your home directory is automatically available.
- The `%test` block in `hapsolo.def` runs smoke tests at the end of the build. If any tool fails to run, the build aborts.
- The image is based on `python:3.11-slim` (~50 MB base) plus the Python dependencies and external binaries. Final image size depends on whether the optional BLAT aligner is included.
