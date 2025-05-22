# Build Docker image
```bash
docker build --no-cache \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  -t seqanalysis:latest .
```

# To start the interactive shell
```bash
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -w /workspace \
  seqanalysis:latest \
  bash
```

# Add the following alias to the .zshrc shell to activate `newenv` from any directory:
```bash
alias seqanalysis='docker run --rm -it -v $(pwd):/workspace -w /workspace seqanalysis:latest'
```

# To run the Makefile
```bash
# Run all FASTQ samples at the same time
make all_samples

# Run a specific FASTQ sample
make CONF=R8FRST_1_384-1x.conf all
make CONF=R8FRST_2_384-2x.conf all
make CONF=R8FRST_3_1536.conf all
make CONF=R8FRST_4_L1-Slow.Birch.conf all
```

# Trim 504 primers from reference sequences (384 & 1536)
```bash
# Trim 504 primers from 384 reference
python scripts/trim_fasta.py refs/lib384_gene_full_with_504.fasta refs/lib384_gene_full_wo_504.fasta

# Trim 504 primers from 1536 reference
python scripts/trim_fasta.py refs/lib1536_gene_full_with_504.fasta refs/lib1536_gene_full_wo_504.fasta
```
