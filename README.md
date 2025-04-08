IN MODIFICATION!
# Y-mer

**Y-mer** is a tool for determining Y chromosome haplogroups from ultra-low-coverage (0.005–1×) sequencing data using Y chromosome-specific k-mers.

Y-mer supports:
- 🧬 Creating and training new models from high-coverage WGS individual reads
- 📦 Using pre-trained models for quick prediction
- 🔍 Accepting either mapped (`.bam`) or unmapped (`.fastq`) data

Results include the most likely haplogroup(s) with statistical support.

---

## 📖 Citation

Puurand T, Möls M, Kaplinski L, Maal K, Krjutskov K, Salumets A, Kivisild T, Remm M. (2025).  
**_Y-mer: A k-mer-based method for determining human Y chromosome haplogroups from ultra-low sequencing depth data_**  
[https://doi.org/10.21203/rs.3.rs-5042960/v1](https://doi.org/10.21203/rs.3.rs-5042960/v1)

---

## 🚀 Quick Start

### Requirements
- `R`
- `perl`
- [`GenomeTester4`](https://github.com/bioinfo-ut/Genometester4) Must be installed and used with the correct path.

### Download Required Files
- Pre-trained models and resources:  
  - [https://doi.org/10.5281/zenodo.15089783](https://doi.org/10.5281/zenodo.15089783)  
  - [https://bioinfo.ut.ee/randomtandem/mudelid/](https://bioinfo.ut.ee/randomtandem/mudelid/)

---

## 🧪 Model Creation & Training

Our pre-trained models are based on European haplogroups (e.g., 1000 Genomes and Estonian Biobank).  
To work with other populations or increase resolution, users can train custom models.

### 1. Prepare Input Files
- Include at least **10 individuals per haplogroup**
- Format for `men.txt` (tab-separated):  
  ```
  HG_ID SAMPLE1_ID SAMPLE2_ID ...
  ```

Modify the existing `men.txt` and `women.txt` files as needed.

### 2. Set Up Directories
```bash
mkdir data lists
```

Edit `Y-mer.pl` to include paths to `.bam` files, then run:
```bash
perl Y-mer.pl
```

### System Requirements
- **SSD**: ~30 GB per sample  
- **RAM**: ~80 GB  
- **Runtime**: ~1.5 hours per sample (SSD speed-dependent)

After processing, retain:
- Model `.Rdata` file  
- `.dbb` database file  
- Final `.txt` count table  
- Optional: male-only k-mer list for future training

> All scripts can be adapted for HPC parallelization (within ~3 hours). We're working on support for this.

---

## 🧬 Haplogroup Prediction (Pre-trained Models)

### Available Models
`M21W`, `M21E`, `M21NE`, `M110W`, `M213E`, `M222NE`, `M43I1`, `M80R1`

Each model uses:
- `.txt` for `glistquery`
- `.dbb` for `gmer_counter`

### 🔗 Model Downloads
- [https://bioinfo.ut.ee/randomtandem/mudelid/](https://bioinfo.ut.ee/randomtandem/mudelid/)  
- [https://doi.org/10.5281/zenodo.15089783](https://doi.org/10.5281/zenodo.15089783)

---

## 🔢 Counting K-mer Frequencies

### Option 1: Using `.fastq` and `gmer_counter`
```bash
gmer_counter -dbb model.dbb /path/sample.fastq | cut -f 3 | tail -n +3 > sample.counts
```

### Option 2: Extract from `.bam`
```bash
samtools fasta sample.bam | gmer_counter -dbb model.dbb - | cut -f 3 | tail -n +3 > sample.counts
```

### Option 3: Using `GenomeTester4` list
```bash
glistquery sample_25.list -f model.txt | cut -f 2 > sample.counts
```

---

## 🔍 Predicting Haplogroups

Run the R script to classify:
```bash
Rscript PREDICTER.R model.Rdata sample.counts sample.Rdata > sample.txt
```

- `sample.txt`: Human-readable haplogroup output  
- `sample.Rdata`: Saved R object for downstream analysis

---

## 🌐 Web Tool

Try the web-based version here:  
🔗 [https://bioinfo.ut.ee/randomtandem/Y-mer/](https://bioinfo.ut.ee/randomtandem/Y-mer/)

- Accepts `.fastq` or `.fastq.gz`  
- Select multiple models  
- Returns `.txt` and `.Rdata` results  
- Max input size: **0.5 GB**

---

## 📁 Example Data

| Sample | Description | Link |
|--------|-------------|------|
| `DA189` | Ancient DNA male sample (Damgaard et al., 2018) | [ERR2505887.fastq.gz](https://bioinfo.ut.ee/randomtandem/mudelid/ERR2505887.fastq.gz) |
| `DA189.bam` | Aligned BAM | [DA189.sort.rmdup.realign.md.bam](https://bioinfo.ut.ee/randomtandem/mudelid/DA189.sort.rmdup.realign.md.bam) |
| `NA20509 chrY` | Assembled chrY (Hallast et al., 2023) | [NA20509.chrY.fasta](https://bioinfo.ut.ee/randomtandem/mudelid/NA20509.HIFIRW.ONTUL.na.chrY.fasta) |

---

## 📬 Contact

For questions or contributions, please open an issue or contact the developers through [bioinfo.ut.ee](https://bioinfo.ut.ee).
