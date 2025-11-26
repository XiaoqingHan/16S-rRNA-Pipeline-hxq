# 16S rRNA Microbiome Analysis Pipeline

A comprehensive pipeline for processing 16S rRNA sequencing data, from raw FASTQ files to taxonomic abundance analysis and statistical comparisons.

## Overview

This pipeline performs end-to-end analysis of paired-end 16S rRNA sequencing data using k-mer based taxonomic classification. It includes quality control, primer trimming, read merging, chimera detection, taxonomic assignment, and statistical analysis.

## Features

- **Quality Control & Preprocessing**: Primer trimming, quality filtering, read merging
- **K-mer Database Construction**: Build custom k-mer databases from reference sequences
- **Optimal K-mer Selection**: Automated selection of optimal k-mer size using Jaccard similarity
- **Chimera Detection**: Remove chimeric sequences using k-mer based approach
- **Taxonomic Classification**: Fast k-mer based taxonomic assignment
- **Statistical Analysis**: Non-parametric testing with FDR correction
- **Visualization**: Composition plots, heatmaps, and boxplots

## Requirements

```
Python >= 3.7
pandas 2.3.0
numpy 1.26.4
scipy 1.11.2
matplotlib 3.9.0
seaborn 0.13.2
statsmodels 0.14.5
sqlite3 3.46.0
```

## Installation

```bash
git clone <repository>
cd microbiome-pipeline
pip install -r requirements.txt
```

## Quick Start

1. **Configure the pipeline** by editing `config.py`:
   ```python
   project = "PRJNA555320"  # Your project name
   groups = ["CF", "SECC"]  # Sample groups
   seq_mode = "paired"      # paired or single
   ```

2. **Run the complete pipeline**:
   ```bash
   # Step 1: Quality control and read merging
   python step1_qc_merge.py
   
   # Step 2: K-mer selection
   python step2_k_selection.py
   
   # Step 3: Chimera detection
   python step3_chimera.py
   
   # Step 4: Taxonomic abundance calculation
   python step4_abundance.py
   
   # Step 5: Statistical analysis and visualization
   python step5_analysis.py
   ```

## Pipeline Steps

### Step 1: Quality Control and Read Merging (`step1_qc_merge.py`)

- Removes primers from 5' end of reads
- Filters reads based on length, quality score, and N content
- Merges paired-end reads with overlap detection
- Outputs: `merged.fasta`, `filtered_R1.fasta`, `filtered_R2.fasta`

**Key Parameters:**
```python
primer_1 = "GTGCCAGCMGCCGCGGTAA"    # Forward primer (V4 region)
primer_2 = "GGACTACHVGGGTWTCTAAT"   # Reverse primer (V4 region)
min_len = 100                       # Minimum sequence length
min_q = 30                          # Minimum quality score
min_overlap = 20                    # Minimum overlap for merging
```

### Step 2: K-mer Selection (`step2_k_selection.py`)

- Tests multiple k-mer sizes (15, 21, 25, 31, 35)
- Calculates Jaccard similarity of top taxa across k values
- Selects optimal k based on consistency and runtime
- Outputs: `kmer_results/kmer_summary.json`

### Step 3: Chimera Detection (`step3_chimera.py`)

- Identifies chimeric sequences using k-mer splitting approach
- Compares left and right halves against reference database
- Removes sequences with conflicting taxonomic assignments
- Outputs: `non_chimeras.fasta`

### Step 4: Taxonomic Abundance (`step4_abundance.py`)

- Assigns taxonomy using optimal k-mer size
- Calculates absolute and relative abundances
- Outputs: `taxon_abundance.csv`

### Step 5: Statistical Analysis (`step5_analysis.py`)

- Generates abundance matrices at genus levels
- Performs non-parametric testing (Mann-Whitney U or Kruskal-Wallis)
- Applies FDR correction for multiple testing
- Creates visualizations (heatmaps, boxplots, composition plots)

## Configuration

### Main Configuration (`config.py`)

```python
# Project settings
root_dir = Path("/path/to/your/data")
project = "PRJNA555320"
groups = ["CF", "SECC"]
seq_mode = "paired"  # or "single"

# Reference database
kmer_db = "/path/to/ref.kmer"
ref_fasta = "/path/to/HOMD_16S_rRNA_RefSeq.fasta"

# Quality control parameters
step1_config = dict(
    primer_1="GTGCCAGCMGCCGCGGTAA",
    primer_2="GGACTACHVGGGTWTCTAAT",
    min_len=100,
    min_q=30,
    max_merge_rate=0.05,
    n_procs=6
)
```

### Database Construction

Build k-mer database from reference sequences:

```bash
python kmer_db.py
```

## Input Data Structure

```
project_directory/
├── Group1/
│   ├── Sample1/
│   │   ├── sample1_1.fastq.gz
│   │   └── sample1_2.fastq.gz
│   └── Sample2/
│       ├── sample2_1.fastq.gz
│       └── sample2_2.fastq.gz
└── Group2/
    └── ...
```

## Output Files

```
project_directory/
├── Group1/
│   └── Sample1/
│       ├── merged.fasta
│       ├── non_chimeras.fasta
│       └── taxon_abundance.csv
├── kmer_results/
│   └── kmer_summary.json
└── abund_results/
    ├── Group1_rel_abund_gen.csv
    ├── heatmap.png
    ├── boxplot.png
    └── non_parametric.csv
```

## Key Features

### K-mer Based Classification
- Fast taxonomic assignment using pre-built k-mer databases
- Automatic optimal k-mer size selection
- Handles IUPAC ambiguous nucleotides

### Statistical Analysis
- Non-parametric tests for group comparisons
- FDR correction for multiple testing
- Focus on user-specified genera of interest

### Visualization
- Stacked bar plots showing microbial composition
- Heatmaps of core genera abundance
- Boxplots for statistical comparisons

## Performance

- Parallel processing support for all computationally intensive steps
- SQLite database for efficient k-mer lookups
- Batch processing for memory efficiency

## Troubleshooting

### Common Issues

1. **Primer not detected**: Check primer sequences and orientation
2. **Low merge rate**: Adjust `max_merge_rate` and `min_overlap` parameters
3. **Memory issues**: Reduce `batch_size` in k-mer analysis

### Logging

All steps generate detailed logs with timestamps and progress information. Check log files for debugging:

```bash
tail -f 16s_analysis_*.log
```

## Citation
Han XQ, (2025). 16S-microbiome-pipeline. 
GitHub repository: https://github.com/XiaoqingHan/16S-rRNA-Pipeline-hxq.git

## Contact

For questions and support, please contact xhan723@hotmail.com.
