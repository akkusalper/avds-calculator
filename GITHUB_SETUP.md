# GitHub Repository Setup Guide

## Initial Setup

### 1. Create New Repository on GitHub

1. Go to [GitHub](https://github.com) and sign in
2. Click the "+" icon in the top right, select "New repository"
3. Repository settings:
   - **Name**: `avds-calculator` (or `variant-detection-score`)
   - **Description**: "Allele-Specific Variant Detection Score - A comprehensive quality assessment metric for genomic variant calling"
   - **Visibility**: Public (recommended for open-source) or Private
   - **DO NOT** initialize with README, .gitignore, or license (we already have these)

### 2. Initialize Local Git Repository

```bash
cd /c/Users/alpyk/OneDrive/Desktop/vLoD_2026

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: Allele-Specific Variant Detection Score Calculator

- Complete mathematical formulation implementation
- Five-component scoring system (VAF*, Q_alt, D*, SB_strict, PB)
- Weighted geometric mean with strand bias artifact detection
- BAM/CRAM and VCF processing pipeline
- Parallel processing support
- VCF annotation output with INFO fields
- Comprehensive testing with synthetic and realistic data
- Docker containerization
- Documentation and examples"
```

### 3. Connect to GitHub and Push

```bash
# Add remote repository (replace [USERNAME] and [REPO-NAME] with your values)
git remote add origin https://github.com/[USERNAME]/[REPO-NAME].git

# Rename branch to main
git branch -M main

# Push to GitHub
git push -u origin main
```

### 4. Repository Configuration (Recommended)

After pushing, configure your repository on GitHub:

#### Topics
Add these topics to help others find your repository:
- `bioinformatics`
- `genomics`
- `variant-calling`
- `quality-control`
- `ngs`
- `sequencing`
- `cancer-genomics`
- `somatic-variants`
- `python`
- `docker`

#### Settings
- Enable **Issues** for bug reports and feature requests
- Enable **Discussions** for Q&A and community support
- Consider enabling **Sponsorships** if you want to accept funding

#### About Section
Add a description in the "About" section:
```
Allele-Specific Variant Detection Score (AVDS) - A training-free, interpretable quality metric for genomic variant calling that combines VAF, quality, depth, strand bias, and position bias into a single score.
```

## Docker Hub Setup (Optional)

### 1. Create Docker Hub Repository

1. Go to [Docker Hub](https://hub.docker.com)
2. Create a new repository named `avds-calculator`
3. Set visibility to Public

### 2. Build and Push Docker Image

```bash
# Login to Docker Hub
docker login

# Build image with your username
docker build -t [YOUR_USERNAME]/avds-calculator:latest .
docker build -t [YOUR_USERNAME]/avds-calculator:v1.0 .

# Push to Docker Hub
docker push [YOUR_USERNAME]/avds-calculator:latest
docker push [YOUR_USERNAME]/avds-calculator:v1.0
```

### 3. Update README

After publishing to Docker Hub, update the README.md Docker section:

Replace:
```bash
docker pull [username]/avds-calculator:latest
```

With:
```bash
docker pull [YOUR_USERNAME]/avds-calculator:latest
```

## Verification Checklist

Before pushing, verify:

- [ ] README.md title is "Allele-Specific Variant Detection Score Calculator"
- [ ] LICENSE file has your name (replace `[Your Name]`)
- [ ] Dockerfile MAINTAINER has your information (replace `[Your Name/Email]`)
- [ ] .gitignore excludes large data files (*.bam, *.vcf.gz, etc.)
- [ ] Test data works: `python test_data/create_test_data.py`
- [ ] Pipeline works: `python avds_pipeline.py -v test_data/test.vcf.gz -b test_data/test.bam -o test_output.tsv --no-parallel`
- [ ] Docker builds: `docker build -t avds-calculator .`
- [ ] Docker runs: `docker run avds-calculator --help`
- [ ] All documentation is in English
- [ ] No Turkish text in any files

## Post-Push Tasks

After successfully pushing to GitHub:

1. **Add a release**: Create v1.0 release with release notes
2. **Update citation**: Add DOI if you create a Zenodo archive
3. **Documentation**: Consider adding a Wiki for detailed documentation
4. **CI/CD**: Set up GitHub Actions for automated testing (optional)
5. **Badge**: Add GitHub stars badge to README
6. **Community**: Add CODE_OF_CONDUCT.md if you expect contributions

## Example GitHub URLs

Replace with your actual URLs:
- Repository: `https://github.com/[USERNAME]/avds-calculator`
- Issues: `https://github.com/[USERNAME]/avds-calculator/issues`
- Releases: `https://github.com/[USERNAME]/avds-calculator/releases`

## Common Issues

### Large File Warning
If you accidentally committed large BAM/VCF files:
```bash
# Remove from git history (use carefully!)
git filter-branch --tree-filter 'rm -f path/to/large/file' HEAD
```

### Authentication Issues
If push fails with authentication error:
```bash
# Use personal access token instead of password
# Generate token at: https://github.com/settings/tokens
# Use token as password when prompted
```

### Line Ending Issues (Windows)
```bash
# Configure git to handle line endings
git config --global core.autocrlf true
```

## Support

For issues with the tool itself, open an issue on GitHub.
For Git/GitHub questions, see [GitHub Docs](https://docs.github.com).
