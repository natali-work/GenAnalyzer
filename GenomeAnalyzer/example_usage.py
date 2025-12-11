#!/usr/bin/env python3
"""
Example usage of the Genome Uncharacterized Gene Analyzer

This script demonstrates how to use the analyzer programmatically.
"""

from genome_uncharacterized_analyzer import UncharacterizedGeneAnalyzer
from Bio import SeqIO

def example_analysis():
    """Example of analyzing a genome."""
    
    # Initialize analyzer with your email
    analyzer = UncharacterizedGeneAnalyzer(
        email="your.email@example.com",  # Replace with your email
        api_key=None  # Optional: add NCBI API key for higher rate limits
    )
    
    # Load genome file
    genome_file = "example_genome.fasta"  # Replace with your genome file
    
    try:
        # Load genome
        records = analyzer.load_genome(genome_file, file_format="fasta")
        
        # Optionally load annotations
        annotations = analyzer.load_annotations("annotations.gff")  # Optional
        
        # Extract uncharacterized genes
        genes = analyzer.extract_genes_from_genome(records, annotations)
        
        print(f"Found {len(genes)} uncharacterized genes")
        
        # Analyze first few genes
        for i, gene in enumerate(genes[:5]):  # Analyze first 5 genes
            print(f"\nAnalyzing gene {i+1}/{min(5, len(genes))}...")
            analysis = analyzer.analyze_gene(gene)
            
            if analysis:
                print(f"Gene: {analysis['gene_id']}")
                print(f"Predicted Function: {analysis['predicted_function']}")
                if analysis['blast_hits']:
                    print(f"Best BLAST hit: {analysis['blast_hits'][0]['description']}")
                    print(f"Identity: {analysis['blast_hits'][0]['identity']:.2f}%")
        
        # Generate full report
        if genes:
            analyzer.generate_report(
                [analyzer.analyze_gene(g) for g in genes[:10]],  # First 10 genes
                "example_report.txt"
            )
    
    except FileNotFoundError:
        print(f"Error: Genome file '{genome_file}' not found.")
        print("Please provide a valid genome file path.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    example_analysis()

