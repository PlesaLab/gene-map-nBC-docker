# DropSynth-Assembled Gene Mapping Pipeline

A reproducible, Docker-based workflow for mapping DropSynth-assembled gene sequences to their reference sequences. The pipeline is orchestrated by a Makefile, parameterized by sample-specific `.conf` files, and runs entirely in a reproducible Docker environment.

---

## 🔍 Overview

This pipeline performs:

1. **Setup:** Creates per-sample output folders (`out/<basename>/`) and a matching log folder (`logs/<basename>/`).
2. **FASTQ QC:** Calculates total read count and length distribution (mean, median, min, max) for the raw FASTQ/FASTQ.GZ.
3. **Motif Processing:** Extracts and processes sequences using a defined primer motif (optionally splits input into chunks).
4. **Optional Trimming:** If `TRIM_SEQS := yes` is set in the `.conf`, runs `trim_padding.py` to remove padding between RE sites and writes a `_nPad.fasta`.
5. **Alignment:** Maps both original and trimmed sequences (if present) to the reference genome using:

   * **BBMap** (`bbmap.sh`)
   * **minimap2** (`minimap2`)
6. **Top Read Extraction:** Extracts top-aligned reads from each alignment result for downstream QC and analysis.

---

## 🚀 Installation & Setup

**1. Clone the Repo**

```bash
git clone git@github.com:PlesaLab/gene-map-nBC-docker.git
cd gene-map-nBC-docker
```

**2. Build the Docker Image**

```bash
docker build --no-cache \
  --platform=linux/arm64 \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  -t newenv:arm64 .
```

Key flags:

* `--no-cache`: clean rebuild
* `--platform`: match your architecture (`linux/arm64` or `linux/amd64`)
* `USER_ID`/`GROUP_ID`: map container permissions to host user
* `-t newenv:arm64`: tags the image

---

## 📦 Data Preparation

Place your **FASTQ/FASTQ.GZ** files into `fastq/` and their matching `.conf` configuration files into `config/`.

Each `.conf` file controls one run and specifies:

* Input FASTQ
* Reference genome
* Barcoding parameters
* Whether to split into chunks (`CHUNK_FASTQ := no`)
* Whether to run trimming on sequence padding (`TRIM_SEQS := yes/no`)

---

## 💻 Running the Pipeline

Run everything from inside the container for reproducibility:

```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  newenv:arm64 \
  bash
```

Once inside the container:

### (OPTIONAL) Combine multiple `FASTQ` replicate files from same library:

```bash
# Twist Unique Libraries
cat fastq/Twist_Unique/replicates/SV825S_1_384-unique1.fastq \
    fastq/Twist_Unique/replicates/SV825S_2_384-unique2.fastq > \
    fastq/Twist_Unique/SV825S_384_unique_combined.fastq

cat fastq/Twist_Unique/replicates/SV825S_3_1536-unique1.fastq \
    fastq/Twist_Unique/replicates/SV825S_4_1536-unique2.fastq > \
    fastq/Twist_Unique/SV825S_1536_unique_combined.fastq
```

### Map a Specific Sample to its Reference

```bash
make CONF=<subdir>/<sample>.conf
```

### Examples:

```bash
# Twist 384 Unique-Overlaps Test Library
make CONF=Twist_Unique/SV825S_Twist_Unique_384.conf

# Twist 1536 Unique-Overlaps Test Library
make CONF=Twist_Unique/SV825S_Twist_Unique_1536.conf
```

You can repeat `make CONF=...` for every dataset you want to process.

Each run generates:

* `out/<basename>/...` — per-sample outputs
* `logs/<basename>/...` — all logs for that sample

### (OPTIONAL) Trim Reference FASTA Files

Trim reference file to remove 50X primers:

- `trim_50X_fasta.py` - standalone python script to remove 50X primers from reference `.fasta` files

```bash
# Trim 504 primers from Twist 384 Unique-Overlaps reference .fasta
python scripts/trim_50X_fasta.py refs/lib384_gene_full_with_504.fasta refs/lib384_gene_full_wo_504.fasta

# Trim 504 primers from Twist 1536 Unique-Overlaps reference .fasta
python scripts/trim_50X_fasta.py refs/lib1536_gene_full_with_504.fasta refs/lib1536_gene_full_wo_504.fasta
```

---

## 🗂️ Configuration File Format

Example `SV825S_Twist_Unique_384.conf`:

```make
INPUT_FASTQ   := $(INPUT_DIR)/Twist_Unique/SV825S_384_unique_combined.fastq
REF_GENOME    := $(REF_DIR)/Twist_Unique/unique_overlap_lib384_gene_full_wo_504.fasta
PROTEIN_FASTA := 

startsitelength := 12
barcode_length  := 20
motif           := ((GATGG)(.)(AACTAACG)){e<2}
start_site      := TGGTAACTAACG
end_site        := ACCAACGGACAA

# Toggle FASTQ chunking: yes or no
# --set to 'no' to skip chunking
CHUNK_FASTQ := no

# Toggle trimming of recovered sequences to remove padding (yes/no)
# --set to 'no' to skip trimming
TRIM_SEQS := no
```

Key variables:

* `INPUT_FASTQ` – path to input directory (`fastq/`) with raw sequence files
* `REF_GENOME` – path to reference directory (`refs/`) with reference sequence files
* `PROTEIN_FASTA` – path to reference directory (`refs/`) with reference protein files
    - (default: blank)
* `startsitelength` – number of base pairs in start site sequence 
* `barcode_length` – number of random nucleotides to add as barcodes to raw sequences
* `motif` – specifies the `startsite` sequence, `bc`, `endsite` sequence, and allowable errors
    - ((start_site)(barcode}(end_site)){e<2>}
* `start_site` - specifies the `start_site` nucleotide sequence
* `end_site` - specifies the `end_site` nucleotide sequence
* `CHUNK_FASTQ` - specifies if the raw .fastq.gz should be split into 1GB chunks to extract BCs
	- yes/no (split input into 1GB chunks)
* `TRIM_SEQS` - specifies if the raw .fastq.gz sequences have additional sequencing padding that needs to be trimmed
	- yes/no (enable trimming step)

---

## 🏗️ Makefile Targets

* **prepare** — creates `out/<basename>/` and `logs/<basename>/`
* **fastq\_length\_stats** — calculates read length stats
* **process\_barcodes** — extracts barcodes (handles chunking if enabled)
* **trim\_padding** — runs only if `TRIM_SEQS := yes`
* **run\_bbmap / run\_minimap** — align reads to reference genome
  (if trimming enabled: also align `_nPad.fasta` and produce separate outputs)
* **extract\_top\_align\_reads\_bbmap / minimap** — extract top reads
  (also runs trimmed equivalents if `TRIM_SEQS := yes`)
* **all** — runs the entire pipeline for the given `.conf`

---

## 📊 Logs & Outputs

* **Logs:** `logs/<basename>/` contains one `.log` per step (fastq stats, barcode processing, BBMap, minimap2, extraction).
* **Outputs:** `out/<basename>/` contains FASTAs, SAM files, counts, and extracted read FASTAs.

---

## 📁 Directory Structure

```bash
├── Dockerfile                    ← default docker settings
├── newenv.yml                    ← default environment file
├── Makefile                      ← default Makefile targets
├── config                        ← default configuration files
│ ├──<basename>/
│ ├──── <basename>.conf
├── fastq/                        ← .fastq.gz-formatted raw sequences
│ ├──<basename>/
│ ├──── <basename>.fastq.gz
├── refs/                         ← .fasta-formatted reference sequences
│ ├──<basename>/
│ ├──── <basename>.lib[]_gene_full_wo_504.fasta
├── scripts/                      ← necessary .py scripts
│ ├── 1_sam_read.R
│ ├── barcode_processing.py
│ ├── extract_top_align_reads.py
│ ├── fastq_length_stats.py
│ ├── split_script.py
│ ├── trim_50X.fasta.py
```

---

## 🛠️ Tips

* Modify `.conf` files to add new datasets.
* Set `TRIM_SEQS := yes` for datasets that need padding removal.
* Inspect `logs/<basename>/` for troubleshooting failed runs.
* Rebuild the Docker image when dependencies or scripts change.

**(Optional) Quick-Access Docker Alias**

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
make CONF=Unique_Twist_384-1536/27BX2S_384_unique_deep.conf
```

7. **To exit the container, run:**
```bash
exit
```

---

## ⚙️ Maintainers
 
- Karl Romanowicz (krom@uoregon.edu)
- Calin Plesa (calin@uoregon.edu)
