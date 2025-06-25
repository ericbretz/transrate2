<p align="center">
  <img src="https://camo.githubusercontent.com/87d2f7858eb2f36b36097cf257faaa40d3146b708d0afb3ebca17f0f2f99ab11/68747470733a2f2f692e696d6775722e636f6d2f31744a76696f442e706e67" alt="TransRate2" width="1000">
</p>

<p align="right">EC Bretz</p>

**TransRate2** is a tool for quality analysis of de-novo transcriptome assemblies. It evaluates assemblies using a wide range of metrics, leveraging read mapping, reference-based statistics, and contig-level features to provide a detailed assessment of transcriptome assembly quality. TransRate2 supports both single and paired-end reads, multiple aligners, and reference-based evaluation.

**TransRate2** is a reimplementation of the original bioinformatics tool, TransRate (Smith-Unna et al. 2015).
https://github.com/blahah/transrate


> [!CAUTION]
> $\Huge\textcolor[RGB]{248, 82, 73}{\textsf{TransRate2 is currently under development}}$

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Dependencies</h2>

- **Python 3.8+**
- **pandas**, **numpy**, **pysam**, **psutil** (Python libraries)
- **Bowtie2** or **Hisat2** - Read aligners
- **Salmon** - Quantification
- **Samtools** - BAM/SAM processing
- **Diamond** (optional, for reference-based analysis)

#### Ubuntu/Debian:
```bash
sudo apt-get install bowtie2 hisat2 salmon samtools diamond
pip install pandas numpy pysam psutil
```

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Installation</h2>

1. Clone the repository:
```bash
git clone https://github.com/ericbretz/transrate2.git
cd TransRate2
```

2. (Optional) Set up a virtual environment and install Python dependencies:
```bash
python3 -m venv venv
source venv/bin/activate
pip install pandas numpy pysam psutil
```

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Usage</h2>

### Basic Usage

```bash
python transrate2.py --assembly assembly.fasta --left reads_1.fq --right reads_2.fq
```

### Command Line Options

| Option         | Short | Description                                      | Default         |
|----------------|-------|--------------------------------------------------|-----------------|
| `--assembly`   | `-a`  | Assembly file(s), comma separated                | *Required*      |
| `--left`       | `-l`  | Left reads file (FASTQ)                          |                 |
| `--right`      | `-r`  | Right reads file (FASTQ)                         |                 |
| `--reference`  | `-f`  | Reference file (FASTA, optional)                 |                 |
| `--output`     | `-o`  | Output directory                                 | `transrate2`    |
| `--threads`    | `-t`  | Number of threads                                | `1`             |
| `--hisat2`     | `-s`  | Use Hisat2 aligner (default: Bowtie2)            | `False`         |
| `--bowtie2`    | `-b`  | Use Bowtie2 aligner                              | `True`          |
| `--clutter`    | `-c`  | Remove intermediate files                        | `False`         |
| `--quiet`      | `-q`  | Suppress terminal output                         | `False`         |
| `--help`       | `-h`  | Display help message                             |                 |

> **Note:** You must provide at least an assembly file. For paired-end reads, provide both `--left` and `--right`. For single-end, provide only one. Reference-based analysis is optional.

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Examples</h2>

#### Basic run with paired-end reads:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq
```

#### Single-end reads:
```bash
python transrate2.py -a assembly.fasta -l reads.fq
```

#### Specify output directory and use Hisat2:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -o results -s
```

#### Reference-based evaluation:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -f reference.fasta
```

#### Remove clutter from output:
```bash
python transrate2.py -a assembly.fasta -l reads_1.fq -r reads_2.fq -c
```

#### Show help:
```bash
python transrate2.py -h
```

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Overview</h2>

1. **Read Mapping** - Maps reads to the assembly using Bowtie2 or Hisat2.
2. **Quantification** - Uses Salmon for transcript quantification.
3. **Contig Metrics** - Calculates a wide range of contig-level and assembly-level metrics (length, GC, N50, coverage, mapping, fragmentation, segmentation, etc.).
4. **Reference-based Evaluation** (optional) - Uses Diamond to compare assembly to a reference and computes reciprocal best hits and coverage.
5. **Scoring** - Computes overall and optimal assembly scores.
6. **Reporting** - Outputs results and summary statistics to the specified directory.

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Output</h2>

TransRate2 creates an output directory (default: `transrate2/`) containing:

```
transrate2/
├── transrate2_<assembly_name>/      # Main results for each assembly
│   ├── <assembly_name>.transrate2.csv   # Main summary metrics
│   ├── <assembly_name>.contigs.csv      # Per-contig metrics
│   ├── good.<assembly_name>.fa          # Good contigs (FASTA)
│   ├── bad.<assembly_name>.fa           # Bad contigs (FASTA)
│   ├── assembly_score_optimisation.csv  # Score optimization data
│   ├── reference/                      # Reference-based results (if used)
│   │   ├── reciprocal_hits.csv
│   │   └── ...
│   └── ...
├── <aligner>_<assembly_name>/       # Aligner-specific files (indices, BAM/SAM)
├── salmon_<assembly_name>/          # Salmon quantification files
├── logs_<assembly_name>/            # Log files
└── assembly.csv                     # Combined summary for all assemblies
```

<h2><img src="https://i.imgur.com/3UA4xwp.png" width="20" align="top">&ensp;Output Files</h2>

- **`<assembly_name>.transrate2.csv`** - Main summary metrics for the assembly
- **`<assembly_name>.contigs.csv`** - Per-contig metrics
- **`good.<assembly_name>.fa`** - Good contigs (FASTA)
- **`bad.<assembly_name>.fa`** - Bad contigs (FASTA)
- **`assembly_score_optimisation.csv`** - Score optimization data
- **`reciprocal_hits.csv`** - Reference-based reciprocal best hits (if reference is provided)
- **`assembly.csv`** - Combined summary for all assemblies
---

**TransRate2** provides a robust, multi-metric assessment of transcriptome assemblies, supporting both reference-free and reference-based evaluation, and is suitable for a wide range of transcriptomic projects.
