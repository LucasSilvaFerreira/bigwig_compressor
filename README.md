# bigwig_compressor

Utilities for attenuating high-expression genes in an RNA bigWig track while keeping the result as a bigWig.

## Install

Clone the repository, then install it with `pip` from the repo root:

```bash
python -m pip install .
```

This installs the CLI command:

```bash
attenuate-rna-bigwig --help
```

### Install directly from GitHub

```bash
python -m pip install "git+https://github.com/LucasSilvaFerreira/bigwig_compressor.git#egg=bigwig-compressor"
```

If you only want the dependencies without installing the package entry point, use:

```bash
python -m pip install -r requirements.txt
```

## Usage

```bash
attenuate-rna-bigwig \
  --input-bigwig input.bw \
  --gtf genes.gtf.gz \
  --output-bigwig output.attenuated.bw \
  --report-dir reports \
  --threshold 50
```

## Outputs

The script writes:

- An attenuated bigWig file
- `gene_signal_report.tsv`
- `threshold_summary.tsv`
- `outlier_threshold_plot.png`
- `attenuated_gene_examples/`
- `run.log`

## Development

For editable installs during local development:

```bash
python -m pip install -e .
```
