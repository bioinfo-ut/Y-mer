IN MODIFICATION!
# Y-mer

**Y-mer** is a tool for determining Y chromosome haplogroups from ultra-low-coverage (0.005–1×) sequencing data using Y chromosome-specific k-mers.

Y-mer supports:
- 🧬 Creating and training new models from high-coverage WGS individual reads
- 📦 Using pre-trained models for quick prediction
- 🔍 Accepting either mapped (`.bam`) or unmapped (`.fastq`) data

Results include the most likely haplogroup(s) with statistics.

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
- [`GenomeTester4`](https://github.com/bioinfo-ut/Genometester4) Must be installed and used with the local installation path.

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

### Available Models .Rdata
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
gmer_counter -dbb model.dbb sample.fastq | cut -f 3 | tail -n +3 > sample.counts
```

### Option 2: Extract from `.bam`
```bash
samtools fasta sample.bam | gmer_counter -dbb model.dbb - | cut -f 3 | tail -n +3 > sample.counts
```

### Option 3: Using `GenomeTester4` listmaker and glistquery
```bash
glistmaker sample.fastq -w 25 -o sample
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

Y-mer uses Y chromosome-specific k-mers and distance-based models to predict Y chromosome haplogroups (Yhg). With this tool the user can upload their own data in the form of a fastq file.  Y-mer will determine the closest Yhg for the uploaded sample in the chosen model on the basis of highest similarity.

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

The distance between sample and haplogroup hg profile, dh, indicates how similar the number of repeats profile of a sample are to the haplogroup average copy-number profile. Look closer output description in supplement methods file for Y-mer publication (Puurand et al. 2024).

The currently available models include 11 basic haplogroups (AB, C, E, G, H, IJ, LT, N, O, Q, R) that are common  at the World (W), 22 (AB, C, E1, E2, E4, G, H, I1, I2, J1, J2, LT, N3, N4, O1, O2'5, O3, O6, Q, R1a, R1b, R2) at European (E), and 23 (E2a, G2a, I1a, I1d, I1i, I1m, I2, Ic, J1, J2a, J2b, LT, N3a3, N3a4, Q, R1a1, R1a2, R1b1, R1b11, R1b2, R1b3, R1b6, R1b8) at Northeast European (NE) levels. The k-mers used in the models have been extracted from sets of 21, 110, 213 and 222 Y chromosomes and the models have been trained on subsets of individuals from the 1000G and EGC projects data. The I1 and R1 models predict only the specified subclades of the given haplogroups.

```bash
[1] "Diagnostics plots will not be produced"
[1] "No sample ID-s found"
[1] "Sample coverages:"
        V1
0.00553951
[1] "Applying CG-related corrections"
[1] "Sequencing coverage uniformity, Dli (smaller is better):"
[1] 0
[1] "Raw distances:"
         AB        C        E       G       H       IJ       LT        N
V1 3.246013 3.214565 3.166282 3.10242 3.39008 3.317615 2.975294 5.084806
         O       Q          R
V1 3.40771 3.08243 -0.4299814
[1] "Most likely haplogroups (assuming contemporary DNA):"
   sample   coverage haplogroup       pvalue alternatives
DA189     DA189 0.00553951          R 1.549354e-17             
```

---

## 📬 Contact

For questions or contributions, please open an issue or contact the developers through [bioinfo.ut.ee](https://bioinfo.ut.ee).
Formatted by ChatGPT
