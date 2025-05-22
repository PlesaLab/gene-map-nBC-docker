# DropSynth-Assembled Gene Mapping Pipeline

## PART I: Setting Up Docker Image for Reproducible Environment

### Prerequisites

- **Git** (to clone the repo)  
- **Docker** (to build and run the container)

---

### 1. Clone the repository

```bash
git clone git@github.com:SynPlexity/gene-map-nBC-docker.git
cd gene-map-nBC-docker
```

If you prefer HTTPS, use:

```bash
git clone https://github.com/SynPlexity/gene-map-nBC-docker.git
cd gene-map-nBC-docker
```

---

### 2. Build the Docker Image

> [!NOTE]
> Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is running in the background (and set default memory to max)

```bash
docker build --no-cache --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) --platform=linux/amd64 -t seqanalysis:latest .
```

- `--no-cache` forces a clean build
- `--build-arg USER_ID` Adds your user id to the image and updates mambauser to make changes to your mounted volume
- `--build-arg GROUP_ID` Same as above, but adds mambauser to the same group as well
- `-t seqanalysis:latest` tages the image for easy reference

---

### 3. Run the Pipeline

All code lives under `/workspace` inside the container.

**Option 1:**

To drop into an interactive shell with your `dropsynth-env:latest` image (and have your `dropsynth` conda env auto‑activated), run:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  seqanalysis:latest \
  bash
```

- `--rm`: remove the container when you exit

- `-it`: interactive TTY

- `-v "$(pwd)":/workspace`: mount your current directory into `/workspace` in the container

- `-w /workspace`: set the working directory inside the container

- `seqanalysis:latest`: the image you just built

- `bash`: start a shell (your `seqanalysis` Conda env will be auto‑activated on launch)

Once the `seqanalysis:latest` image is activated, run `make` targets as follows:

```bash
# To run a specific FASTQ sample
make CONF=R8FRST_1_384-1x.conf all
make CONF=R8FRST_2_384-2x.conf all
make CONF=R8FRST_3_1536.conf all
make CONF=R8FRST_4_L1-Slow.Birch.conf all

# To run all FASTQ samples at the same time
make all_samples
```

**Option 2:**

To mount your project directory and run your full workflow on a specific FASTQ, run:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  seqanalysis:latest \
  mmake CONF=R8FRST_1_384-1x.conf all
```

**Option 3:**

To mount your project directory and run your full workflow on all FASTQ samples, run:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  seqanalysis:latest \
  make all_samples
```

---

### 4. One-Line Shortcut (Optional)

Clone, build, and run everything in a single command:

```bash
git clone git@github.com:SynPlexity/gene-map-nBC-docker.git \
  && cd gene-map-nBC-docker \
  && docker build -t seqanalysis:latest . \
  && docker run --rm -it -v "$(pwd)":/workspace -w /workspace seqanalysis:latest make all_samples
```

---

### 5. Add Alias for Quick-Access to Container (Optional)

> [!NOTE]
> The following option is intended for Mac users only

Add an alias to your `~/.zshrc`, then reload it:

1. **Edit your** `~/.zshrc:`
```bash
nano ~/.zshrc
```

2. **Add this line** (you can also include --rm so containers auto‑cleanup on exit):
```bash
alias seqanalysis='docker run --platform=linux/amd64 --rm -it -v $(pwd):/workspace -w /workspace seqanalysis:latest'
```

3. **Save & exit** (`Ctrl+O`, then `Enter`, then `Ctrl-X`).

4. **Reload your shell config:**
```bash
source ~/.zshrc
```

5. **Now, anywhere inside your project directory you can just run:**
```bash
seqanalysis
```

6. **With the container running in the `/workspace` directory, run `Make` targets:**
```bash
make all_samples
```

7. **To exit the container, run:**
```bash
exit
```

---

## PART II: DropSynth Oligo-Generation Pipeline

> [!Warning]
> This pipeline is still in beta - some features have not been exhaustively tested.

### Overview

This repository packages a Makefile‑driven pipeline (inside a Docker container) for mapping DropSynth-assembled gene sequencing reads (`fastq.gz`) to their reference gene sequences (`fasta`).

It does the following:

1. Initializes directory structure for logs and outputs
2. Extracts and processes barcode sequences from FASTQ input
3. Aligns sequences to reference genome using BBMap or MiniMap2
4. (Optional) Extracts top-aligned reads for downstream analysis

Everything is orchestrated via a single `Makefile`, configured by individual `*.conf` files specific for each assembled library (`fastq.gz`), and runnable in a reproducible Docker environment.

---

### Repository Layout

```bash
├── Dockerfile                    ← default docker settings
├── seqanalysis.yml               ← default environment file
├── Makefile                      ← default Makefile targets
├── config                        ← default configuration files
│ ├── [fastq.gz.specific_1].conf
│ ├── [fastq.gz.specific_2].conf
│ ├── [fastq.gz.specific_3].conf
├── fastq/                        ← .fastq.gz-formatted raw sequences
│ ├── [fastq.gz.specific_1].fastq.gz
│ ├── [fastq.gz.specific_2].fastq.gz
│ ├── [fastq.gz.specific_3].fastq.gz
├── refs/                         ← .fasta-formatted reference sequences
│ ├── lib[]_gene_full_wo_504.fasta
│ ├── custom[]_out_split_Lib[].full_nRE_nPrim.genes
├── scripts/                      ← necessary .py scripts
│ ├── trim_fasta.py
│ ├── barcode_processing.py
│ ├── extract_top_align_reads.py
```

---

### Configuration (`[fastq.gz.specific_1].conf`)

All parameters live in a single `[fastq.gz.specific_1].conf`.  Key variables include:

- `INPUT_FASTQ` – path to input directory (`fastq/`) with raw sequence files
- `REF_GENOME` – path to reference directory (`refs/`) with reference sequence files
- `PROTEIN_FASTA` – path to reference directory (`refs/`) with reference protein files
    - (default: blank)
- `startsitelength` – number of base pairs in start site sequence 
- `barcode_length` – number of random nucleotides to add as barcodes to raw sequences
- `motif` – specifies the `startsite` sequence, `bc`, `endsite` sequence, and allowable errors
    - ((start_site)(barcode}(end_site)){e<2>}
- `start_site` - specifies the `start_site` nucleotide sequence
- `end_site` - specifies the `end_site` nucleotide sequence

You can duplicate `[fastq.gz.specific_1].conf` for multiple libraries.

---

### Makefile Targets (`Makefile`)

- `prepare` – create log and output directories

- `process_barcodes` – `barcode_processing.py` → Extract BCs from input fastq.gz file

- `run_bbmap` – `bbmap.sh` → Align and map extracted sequences to reference sequences (BBMap)

- `run_minimap` – `minimap2` → Align and map extracted sequences to reference sequences (MiniMap2)

- `extract_top_align_reads_bbmap` – `extract_top_align_reads.py` → Extract top aligned reads from BBMap

- `extract_top_align_reads_minimap` – `extract_top_align_reads.py` → Extract top aligned reads from MiniMap2

- `all` – run all targets in sequential order for a specific `[fastq.gz.specific_1].conf` dataset

- `all_samples` – run all targets in sequential order for all `[fastq.gz.specific_X].conf` datasets

**OPTIONAL**

- `trim_fasta.py` - standalone python script to remove 504 primers from reference `.fasta` files

```bash
# Trim 504 primers from 384 reference
python scripts/trim_fasta.py refs/lib384_gene_full_with_504.fasta refs/lib384_gene_full_wo_504.fasta

# Trim 504 primers from 1536 reference
python scripts/trim_fasta.py refs/lib1536_gene_full_with_504.fasta refs/lib1536_gene_full_wo_504.fasta
```

---

### Logs & Outputs

- **Logs:** each script writes to `logs/<step>.log`
- **Results:** look under `out/` for subdirectories by step

---

### Extending & Troubleshooting

- Modify or clone `*.conf` file for new libraries
- Edit `Makefile` if you add new scripts or targets
- Rebuild Docker image (`Dockerfile`) if you change dependencies
- Inspect `logs/` for errors or exceptions

---

**Maintainers:**
 
- Karl Romanowicz (karl@synplexity.com)
- Calin Plesa (calin@synplexity.com)
