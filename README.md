# Allele-Specific Variant Detection Score Calculator

A comprehensive quality assessment metric for genomic variant calling.

## Overview

The Allele-Specific Variant Detection Score (AVDS) is a training-free, interpretable quality metric that combines five key sequencing quality indicators into a single score (0-100):

- **VAF*** - Normalized Variant Allele Frequency (weight: 0.30)
- **Q̄_alt** - Average Alternative Allele Quality (weight: 0.35)
- **D*** - Normalized Sequencing Depth (weight: 0.10)
- **SB_strict** - Strict Strand Bias with 25% minimum rule (weight: 0.15)
- **PB** - Position Bias (weight: 0.10)

## Installation

```bash
# Clone or download the repository
cd vLoD_2026

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### Basic Usage

```bash
python avds_pipeline.py \
  -v variants.vcf \
  -b alignment.bam \
  -o avds_results.tsv
```

### CRAM Files (requires reference)

```bash
python avds_pipeline.py \
  -v variants.vcf \
  -b alignment.cram \
  -r reference.fa \
  -o avds_results.tsv
```

### Low-Frequency Variants (e.g., ctDNA)

```bash
python avds_pipeline.py \
  -v variants.vcf \
  -b alignment.bam \
  -o avds_results.tsv \
  --theta-vaf 0.10
```

### Parallel Processing

```bash
# Use 8 threads
python avds_pipeline.py \
  -v variants.vcf \
  -b alignment.bam \
  -o avds_results.tsv \
  -t 8
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-v, --vcf` | Input VCF file | Required |
| `-b, --bam` | Input BAM/CRAM file | Required |
| `-r, --reference` | Reference FASTA (for CRAM) | None |
| `-o, --output` | Output TSV file | Required |
| `-s, --sample` | Sample name | First sample |
| `-t, --threads` | Number of threads | CPU count - 1 |
| `--theta-vaf` | VAF normalization threshold | 0.50 |
| `--theta-depth` | Depth normalization threshold | 100 |
| `--min-coverage` | Minimum coverage | 5 |
| `--no-parallel` | Disable parallel processing | False |
| `--verbose` | Enable verbose logging | False |

## Score Interpretation

| Score Range | Quality | Recommended Action |
|-------------|---------|-------------------|
| 90-100 | Excellent | Clinical reporting |
| 70-89 | High Confidence | Research/Clinical |
| 50-69 | Moderate | Research use |
| 30-49 | Low | Manual review required |
| <30 | Artifact | Filter out |

## Output Format

The output TSV file contains:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `pos` | Position (1-based) |
| `id` | Variant ID |
| `ref` | Reference allele |
| `alt` | Alternative allele |
| `avds_score` | Final AVDS score (0-100) |
| `vaf_raw` | Raw VAF |
| `vaf_norm` | Normalized VAF* |
| `q_alt_mean` | Average alt quality |
| `depth_total` | Total depth |
| `depth_norm` | Normalized depth |
| `sb_strict` | Strand bias score |
| `sb_min_ratio` | Minimum strand ratio |
| `sb_failed_25` | Failed 25% rule (True/False) |
| `pb_score` | Position bias score |
| `pb_avg_pos` | Average position (%) |
| `n_alt` | Alternative allele count |
| `n_ref` | Reference allele count |
| `n_forward` | Forward strand count |
| `n_reverse` | Reverse strand count |
| `quality_category` | Quality category |
| `recommended_action` | Recommended action |

## Python API

```python
from avds_pipeline import AVDSPipeline
from avds_calculator import AVDSConfig

# Create configuration
config = AVDSConfig(
    theta_vaf=0.10,  # Low-frequency variant detection
    theta_depth=100,
    min_coverage=10
)

# Initialize pipeline
pipeline = AVDSPipeline(
    bam_path="alignment.bam",
    reference_path="reference.fa",
    config=config,
    n_threads=8
)

# Process VCF
results = pipeline.process_vcf(
    vcf_path="variants.vcf",
    output_path="avds_results.tsv",
    parallel=True
)

# Access results as DataFrame
print(results.head())
print(f"Mean AVDS: {results['avds_score'].mean():.2f}")
```

## Docker Usage

### Pull from Docker Hub (when published)

```bash
docker pull [username]/avds-calculator:latest
```

### Build locally

```bash
cd vLoD_2026
docker build -t avds-calculator .
```

### Run with Docker

```bash
# Basic usage
docker run -v /path/to/data:/data avds-calculator \
  -v /data/variants.vcf.gz \
  -b /data/alignment.bam \
  -o /data/avds_results.tsv

# With custom parameters
docker run -v $(pwd)/data:/data avds-calculator \
  -v /data/variants.vcf.gz \
  -b /data/alignment.bam \
  -o /data/avds_results.tsv \
  --theta-vaf 0.10 \
  -t 8
```

### Using docker-compose

```bash
# Edit docker-compose.yml with your file paths
docker-compose up
```

## Version History

- **v1.0** (February 2026) - Initial release
  - Complete AVDS mathematical formulation
  - Parallel processing support
  - BAM/CRAM compatibility
  - VCF annotation output
  - Docker containerization
