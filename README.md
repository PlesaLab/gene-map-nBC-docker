# DropSynth-Assembled Gene Mapping Pipeline

A reproducible, Docker-based workflow for mapping DropSynth-assembled gene sequences to their reference sequences, orchestrated by a Makefile and configured with sample-specific `.conf` files.

## 🔍 Overview

This workflow does the following:

1. Initializes directory structure for logs and outputs
2. Extracts and processes barcode sequences from FASTQ input
3. Aligns sequences to reference genome using BBMap and MiniMap2
4. Extracts top-aligned reads for downstream analysis

The entire workflow is orchestrated through a single `Makefile`, with individual `.conf` files specifying parameters for each assembled library (`fastq.gz`). Designed for reproducibility, the pipeline runs within a Docker environment and processes all input `fastq.gz` files sequentially using the `make all_samples` target—provided each file has a corresponding `.conf` configuration.

## 🚀 Installation & Setup

**1. Clone the Repo**

```bash
git clone git@github.com:PlesaLab/gene-map-nBC-docker.git
cd gene-map-nBC-docker
```

If you prefer HTTPS, use:

```bash
git clone https://github.com/PlesaLab/gene-map-nBC-docker.git
cd gene-map-nBC-docker
```

**2. Build the Docker Image**

> [!NOTE]
> Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is running in the background (and change default memory to max in Settings)

```bash
docker build --no-cache \
  --platform=linux/arm64 \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  -t newenv:arm64 .
```

- `--no-cache` forces a clean build
- `--platform=linux/arm64` specifies native architecture for build (***change as needed***)
- `--build-arg USER_ID` Adds your user id to the image and updates mambauser to make changes to your mounted volume
- `--build-arg GROUP_ID` Same as above, but adds mambauser to the same group as well
- `-t newenv:arm64` tags the image for easy reference

## 📦 Data Preparation

> [!NOTE]
> To run a reproducible example with this pipeline, clone this repository with the included `SV825S_384_unique_combined.fastq` and `SV825S_1536_unique_combined.fastq` **input files** as well as the `SV825S_384_unique_combined.conf` and `SV825S_1536_unique_combined.conf` **configuration files**. This dataset corresponds to the `Twist-Test-Unique-Overlaps` assembled libraries.

All code lives under `/workspace` inside the container.

**1. Run Interactively**

To drop into an interactive shell with your `newenv:latest` image (and have your `newenv` conda env auto‑activated), run:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  newenv:arm64 \
  bash
```

- `--rm`: remove the container when you exit

- `-it`: interactive TTY

- `-v "$(pwd)":/workspace`: mount your current directory into `/workspace` in the container

- `-w /workspace`: set the working directory inside the container

- `newenv:latest`: the image you just built

- `bash`: start a shell (your `newenv` Conda env will be auto‑activated on launch)

Once the `newenv:latest` image is activated, run `make` targets as follows:

***To run a specific FASTQ sample:***

```bash
# Twist 384-1536 Unique-Overlaps Test Libraries
make CONF=Unique_Twist_384-1536/SV825S_384_unique_combined.conf
make CONF=Unique_Twist_384-1536/SV825S_1536_unique_combined.conf
```

***To run all FASTQ samples at the same time:***

```bash
# Twist 384-1536 Unique-Overlaps Test Libraries
make CONFIG_SUBDIR=Unique_Twist_384-1536 all_samples
```

***If you need to combine multiple `FASTQ` replicate files from the same library:***

```bash
# Twist 384-1536 Unique
cat fastq/Unique_Twist_384-1536/replicates/SV825S_1_384-unique1.fastq \
    fastq/Unique_Twist_384-1536/replicates/SV825S_2_384-unique2.fastq > \
    fastq/Unique_Twist_384-1536/SV825S_384_unique_combined.fastq

cat fastq/Unique_Twist_384-1536/replicates/SV825S_3_1536-unique1.fastq \
    fastq/Unique_Twist_384-1536/replicates/SV825S_4_1536-unique2.fastq > \
    fastq/Unique_Twist_384-1536/SV825S_1536_unique_combined.fastq
```

**2. (Optional) Quick-Access Alias**

Add Alias for Quick-Access to Container

> [!NOTE]
> The following option is intended for Mac users only

Add an alias to your `~/.zshrc`, then reload it:

1. **Edit your** `~/.zshrc:`
```bash
nano ~/.zshrc
```

2. **Add this line** (you can also include --rm so containers auto‑cleanup on exit):
```bash
alias newenv='docker run --rm -it -v "$(pwd)":/workspace -w /workspace newenv:arm64'
```

3. **Save & exit** (`Ctrl+O`, then `Enter`, then `Ctrl-X`).

4. **Reload your shell config:**
```bash
source ~/.zshrc
```

5. **Now, anywhere inside your project directory you can just run:**
```bash
newenv
```

6. **With the container running in the `/workspace` directory, run `Make` targets:**
```bash
make all_samples
```

7. **To exit the container, run:**
```bash
exit
```

## 📁 Directory Structure

```bash
├── Dockerfile                    ← default docker settings
├── newenv.yml                    ← default environment file
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

## 🗂️ Configuration Files

All parameters live in a single `[fastq.gz.specific].conf`.  Key variables include:

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
- `CHUNK_FASTQ` - specifies if the raw .fastq.gz should be split into 1GB chunks to extract BCs

You can duplicate `[fastq.gz.specific].conf` for multiple libraries.

## 🏗️ Makefile Targets

- `prepare` – create log and output directories

- `process_barcodes` – `barcode_processing.py` → Extract BCs from input fastq.gz file

    - (OPTIONS) `split_fastq` - `split_script.py` → Split raw `.fastq.gz` into 1GB chunks
    - (OPTIONS) `process_barcodes_chunks` - `process_barcodes` → Extract BCs from input chunks
    - (OPTIONS) `merge_fasta` → Merges `chunk.fasta` into `.fasta` for downstream mapping

- `run_bbmap` – `bbmap.sh` → Align and map extracted sequences to reference sequences (BBMap)

- `run_minimap` – `minimap2` → Align and map extracted sequences to reference sequences (MiniMap2)

- `extract_top_align_reads_bbmap` – `extract_top_align_reads.py` → Extract top aligned reads from BBMap

- `extract_top_align_reads_minimap` – `extract_top_align_reads.py` → Extract top aligned reads from MiniMap2

- `all` – run all targets in sequential order for a specific `[fastq.gz.specific].conf` dataset

- `all_samples` – run all targets in sequential order for all `[fastq.gz.specific].conf` datasets

**OPTIONAL**

- `trim_fasta.py` - standalone python script to remove 504 primers from reference `.fasta` files

```bash
# Trim 504 primers from 384 reference
python scripts/trim_fasta.py refs/lib384_gene_full_with_504.fasta refs/lib384_gene_full_wo_504.fasta

# Trim 504 primers from 1536 reference
python scripts/trim_fasta.py refs/lib1536_gene_full_with_504.fasta refs/lib1536_gene_full_wo_504.fasta
```

## 📊 Logs & Outputs

- **Logs:** each script writes to `logs/<step>.log`
- **Results:** look under `out/` for subdirectories by step

## 🛠️ Extending & Troubleshooting

- Modify or clone `[fastq.gz.specific].conf` file for new libraries
- Edit `Makefile` if you add new scripts or targets
- Rebuild Docker image (`Dockerfile`) if you change dependencies
- Inspect `logs/` for errors or exceptions

## ⚙️ Maintainers
 
- Karl Romanowicz (krom@uoregon.edu)
- Calin Plesa (calin@uoregon.edu)
