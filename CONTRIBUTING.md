# Contributing to AVDS Calculator

Thank you for your interest in contributing to the Allele-Specific Variant Detection Score (AVDS) Calculator! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:
- A clear, descriptive title
- Steps to reproduce the issue
- Expected vs. actual behavior
- Your environment (OS, Python version, package versions)
- Sample data or code to reproduce (if possible)

### Suggesting Enhancements

Feature requests are welcome! Please open an issue with:
- A clear description of the feature
- Use cases and benefits
- Any relevant references or implementations from other tools

### Pull Requests

1. Fork the repository
2. Create a new branch for your feature (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests if applicable
5. Ensure all tests pass
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to your branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Code Style

- Follow PEP 8 guidelines for Python code
- Use meaningful variable and function names
- Add docstrings to functions and classes
- Comment complex logic
- Keep functions focused and modular

### Testing

- Test your changes with both synthetic and real data
- Ensure backward compatibility
- Add unit tests for new features
- Verify Docker builds and runs correctly

## Development Setup

```bash
# Clone your fork
git clone https://github.com/[YOUR_USERNAME]/avds-calculator.git
cd avds-calculator

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests
python test_data/create_test_data.py
python avds_pipeline.py -v test_data/test.vcf.gz -b test_data/test.bam -o test_output.tsv --no-parallel
```

## Project Structure

```
vLoD_2026/
├── avds_calculator.py     # Core AVDS calculation logic
├── avds_pipeline.py       # VCF/BAM processing pipeline
├── requirements.txt       # Python dependencies
├── README.md              # User documentation
├── test_data/             # Synthetic test data
└── public_somatic_test/   # Realistic test data
```

## Mathematical Formulation

The AVDS score is based on a weighted geometric mean of five components:

```
AVDS = 100 × [VAF*^0.30 × Q̄_alt^0.35 × D*^0.10 × SB_strict^0.15 × PB^0.10]
```

When modifying scoring logic:
- Preserve the mathematical formulation
- Maintain interpretability
- Document any changes to weights or thresholds
- Validate with test cases

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Provide constructive feedback
- Focus on the technical merit of ideas

## Questions?

Feel free to open an issue for questions or reach out to the maintainers.

Thank you for contributing!
