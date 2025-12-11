# Genome Uncharacterized Gene Analyzer

## Overview
This script scans a genome, identifies uncharacterized/hypothetical genes, and predicts their function based on:
- BLAST homology searches against NCBI databases
- Protein domain analysis
- Sequence similarity to known proteins

## Installation

```bash
pip install -r requirements_genome_analyzer.txt
```

## Usage

### Basic Usage
```bash
python genome_uncharacterized_analyzer.py genome.fasta
```

### With Annotation File (GFF/GTF)
```bash
python genome_uncharacterized_analyzer.py genome.fasta -a annotations.gff
```

### Specify Email for NCBI (Required for BLAST)
```bash
python genome_uncharacterized_analyzer.py genome.fasta -e your.email@example.com
```

### Full Options
```bash
python genome_uncharacterized_analyzer.py genome.fasta \
    -a annotations.gff \
    -e your.email@example.com \
    -o report.txt \
    --max-genes 20 \
    --min-length 300
```

## Parameters

- `genome_file`: Path to genome file (FASTA, GenBank, or EMBL format)
- `-a, --annotations`: Path to annotation file (GFF/GTF format)
- `-o, --output`: Output report file (default: uncharacterized_genes_report.txt)
- `-e, --email`: Email for NCBI Entrez (required for BLAST searches)
- `-k, --api-key`: NCBI API key (optional, for higher rate limits)
- `-f, --format`: Genome file format (fasta, genbank, embl)
- `--min-length`: Minimum ORF length in nucleotides (default: 300)
- `--max-genes`: Maximum number of genes to analyze (default: 10)

## Output

The script generates:
1. **Text Report** (`uncharacterized_genes_report.txt`): Human-readable analysis report
2. **JSON Report** (`uncharacterized_genes_report.json`): Machine-readable data

## Features

- Identifies uncharacterized genes by keywords (hypothetical, unknown, putative, etc.)
- Finds ORFs in unannotated sequences
- Translates DNA to protein sequences
- Performs BLAST searches against NCBI databases
- Predicts protein function based on homology
- Generates comprehensive reports

## Notes

- BLAST searches require internet connection and may take time
- NCBI email is required for API access
- For large genomes, consider using `--max-genes` to limit analysis
- Local BLAST+ installation can be used instead of web BLAST for faster results

