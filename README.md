# DropSynth Gene Design Pipeline

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
> Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is running in the background

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
> This software is still in beta - some features have not been exhaustively tested.
> Use at your own risk.

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
├── Dockerfile				← default docker settings
├── seqanalysis.yml 	← default environment file
├── Makefile				  ← default Makefile targets
├── config 			      ← default parameter files
│ ├── [fastq.gz.specific_1].conf
│ ├── [fastq.gz.specific_2].conf
│ ├── [fastq.gz.specific_3].conf
├── fastq/ 				    ← .fastq.gz-formatted raw sequence inputs
│ ├── [fastq.gz.specific_1].fastq.gz
│ ├── [fastq.gz.specific_2].fastq.gz
│ ├── [fastq.gz.specific_3].fastq.gz
├── refs/ 				    ← .fasta or .genes-formatted reference sequences
│ ├── lib[]_gene_full_wo_504.fasta
│ ├── custom[]_out_split_Lib[].full_nRE_nPrim.genes
├── scripts/ 				  ← necessary .py scripts
│ ├── trim_fasta.py
│ ├── barcode_processing.py
│ ├── extract_top_align_reads.py
```

---

### Configuration (`Library.conf`)

All parameters live in a single `Library.conf`.  Key variables include:

- `codon_usage_file` – path to a JSON or CSV of codon frequencies  
- `input_type` – `AA` or `DNA`  
- `input_file` – source FASTA/Excel filename  
- `num_oligos` – target oligos per gene  
- cloning sites – `cloning_fwd`, `cloning_rev`, etc.  
- enzyme rules – `REsites_file` pointing at your CSV of restriction enzymes  
- `div_lib_size`, `number_of_libs` – how to chunk split genes into sub‑libraries  
- buffer settings – `padding_var`, `padding_length`, etc.  
- homopolymer thresholds, primer files, barcode sets, etc.

You can duplicate or override `Library.conf` for multiple libraries.

---








# Trim 504 primers from reference sequences (384 & 1536)
```bash
# Trim 504 primers from 384 reference
python scripts/trim_fasta.py refs/lib384_gene_full_with_504.fasta refs/lib384_gene_full_wo_504.fasta

# Trim 504 primers from 1536 reference
python scripts/trim_fasta.py refs/lib1536_gene_full_with_504.fasta refs/lib1536_gene_full_wo_504.fasta
```
