IN MODIFICATION!
# Y-mer

**Y-mer** is a tool for determining Y chromosome haplogroups from ultra-low-coverage (0.005–1× with high confidence) sequencing data using Y chromosome-specific k-mers.

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

## 🚀 Quick Start In Linux

### Requirements
- `R`
- `perl`
- [`GenomeTester4`](https://github.com/bioinfo-ut/Genometester4) Must be installed and used with the local installation path.
- [`samtools`](https://github.com/samtools/samtools) if starting with BAM/CRAM files 

### Download Required Files
It is recommended to create a directory listing in by model ID (M21W, M21E, M21NE, M110W, M213E, M222NE, M43I1, and M80R1), and download the model.Rdata and model.txt or/and model.dbb files by model ID user plan to use.The model.Rdata file contains information about the k-mer frequencies by haplogroups, model.txt and model.dbb files contain the final list of k-mers used in the model. 
- Pre-trained models and resources:  
  - [https://doi.org/10.5281/zenodo.15089783](https://doi.org/10.5281/zenodo.15089783)  
  - [https://bioinfo.ut.ee/randomtandem/mudelid/](https://bioinfo.ut.ee/randomtandem/mudelid/)

---
### File Types & Data Structures
- sample - in most cases, it is a fastq file prefix (ID), but every temporary file is identified by the ID of the sample
- model - mainly Rdata R formatted file prefix (ID), containing information needed for calling HG represented in the model. 
- lists - mainly temporary binary files containing information about k-mer sequences and frequencies. The current workflow contains different k-mer manipulation options to prepare data for the model.
- tables - collected k-mers with frequencies from male samples to inputs for MWT.R and MODEL.R
- temporary files - files either selecting k-mers via list files or used for the generation of table files.
- result files from model training - model.Rdata, model.dbb and model.txt

---
## 🧪 Model Creation & Training

At the moment, our ready-to-use models are adapted for the detection of the main sub-clades of haplogroups common in present-day Europe and miss many important haplogroups that are either uncommon or frequent outside Europe.  These restrictions were set by our use of  the 1000 Genomes Project and the Estonian Biobank data as references in the models we have generated and tested. When working with data from other world regions or when aiming for higher haplogroup resolution within a region, the users can design their own haplogroup lists and train their own models based on high quality reference data they have access to. 

Our pre-trained models are based on European haplogroups (e.g., 1000 Genomes and Estonian Biobank).  
To work with other populations or increase resolution, users can train custom models.

### 1. Prepare Input Files
The first step of creating a new model involves the generation of a list from bam( cram or fastq) files of high quality genomes representing, ideally with at least 10 individuals per each targeted haplogroup, from the range of haplogroups to be examined. The IDs of each of these bam files should be presented as a list in a table, similar to the example file [`men.txt`]([https://github.com/samtools/samtools](https://github.com/bioinfo-ut/Y-mer/blob/main/model_training/men.txt)) . In this tab-separated file, each line represents one haplogroup to be included. The name of each haplogroup is shown in the first column. Other columns show ID-s of individuals from the given haplogroup. There is no limit set to the number of individuals but 10 individuals is advisable as a minimum.
 
The structure of the women.txt file, containing the ID-s of female WGS data from which femal k-mer lists will be created, is the same as the men.txt but has only just one row, where the entry in the first column should be ‘N’, followed by entries of the IDs of female WGS data to be used. In case of the available models, we have used 15 female high-coverage genomes for building female k-mer lists.

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
- Optional: male-only k-mer binary list for future training

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

The distance between the sample and the haplogroup gh profiles, dh, indicates how similar the number of repeats profile of a sample is to the haplogroup average copy-number profile. A more detailed description of the output can be found in the supplementary methods file of the Y-mer publication (Puurand et al. 2025).

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
