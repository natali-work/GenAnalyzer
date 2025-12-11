#!/usr/bin/env python3
"""
Genome Uncharacterized Gene Analyzer

This script scans a genome, identifies uncharacterized genes, and predicts
their function based on homology searches and domain analysis.
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import json
from datetime import datetime

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez
    from Bio import ExPASy
    from Bio.SwissProt import SwissProt
except ImportError:
    print("ERROR: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

try:
    import requests
except ImportError:
    print("ERROR: requests is required. Install with: pip install requests")


class UncharacterizedGeneAnalyzer:
    """Analyzes uncharacterized genes in a genome."""
    
    # Keywords that indicate uncharacterized/hypothetical proteins
    UNCHARACTERIZED_KEYWORDS = [
        'hypothetical', 'uncharacterized', 'unknown', 'putative',
        'predicted', 'similar to', 'conserved', 'domain-containing',
        'orf', 'open reading frame', 'unnamed', 'unnamed protein'
    ]
    
    def __init__(self, email: str = "your.email@example.com", api_key: Optional[str] = None):
        """
        Initialize the analyzer.
        
        Args:
            email: Email for NCBI Entrez (required for API access)
            api_key: Optional NCBI API key for higher rate limits
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.results = []
        
    def load_genome(self, genome_file: str, file_format: str = "fasta") -> List[SeqRecord]:
        """
        Load genome sequences from a file.
        
        Args:
            genome_file: Path to genome file
            file_format: File format (fasta, genbank, etc.)
            
        Returns:
            List of sequence records
        """
        print(f"Loading genome from {genome_file}...")
        try:
            records = list(SeqIO.parse(genome_file, file_format))
            print(f"Loaded {len(records)} sequences")
            return records
        except Exception as e:
            print(f"Error loading genome: {e}")
            sys.exit(1)
    
    def load_annotations(self, annotation_file: Optional[str] = None, 
                        format: str = "gff") -> Dict:
        """
        Load gene annotations from GFF/GTF file.
        
        Args:
            annotation_file: Path to annotation file
            format: File format (gff, gtf)
            
        Returns:
            Dictionary mapping gene IDs to annotation data
        """
        annotations = {}
        if not annotation_file or not os.path.exists(annotation_file):
            print("No annotation file provided. Will extract from sequence headers.")
            return annotations
            
        print(f"Loading annotations from {annotation_file}...")
        try:
            with open(annotation_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        attributes = parts[8]
                        # Extract gene ID and description
                        gene_id = None
                        description = ""
                        for attr in attributes.split(';'):
                            if 'ID=' in attr or 'gene_id=' in attr:
                                gene_id = attr.split('=')[1].strip()
                            if 'Name=' in attr or 'gene=' in attr:
                                description = attr.split('=')[1].strip()
                        if gene_id:
                            annotations[gene_id] = {
                                'description': description,
                                'type': parts[2],
                                'start': int(parts[3]),
                                'end': int(parts[4]),
                                'strand': parts[6]
                            }
            print(f"Loaded {len(annotations)} annotations")
        except Exception as e:
            print(f"Warning: Could not load annotations: {e}")
            
        return annotations
    
    def is_uncharacterized(self, description: str) -> bool:
        """
        Check if a gene/protein is uncharacterized based on its description.
        
        Args:
            description: Gene/protein description
            
        Returns:
            True if uncharacterized
        """
        desc_lower = description.lower()
        return any(keyword in desc_lower for keyword in self.UNCHARACTERIZED_KEYWORDS)
    
    def extract_genes_from_genome(self, records: List[SeqRecord], 
                                  annotations: Dict = None) -> List[Dict]:
        """
        Extract gene sequences from genome.
        
        Args:
            records: Genome sequence records
            annotations: Optional annotation dictionary
            
        Returns:
            List of gene dictionaries
        """
        genes = []
        
        for record in records:
            # Try to extract genes from sequence header
            header = record.description
            gene_id = record.id
            description = header
            
            # Check if uncharacterized
            if self.is_uncharacterized(description):
                # Try to find ORFs in the sequence
                seq = str(record.seq)
                
                # Simple ORF finding (minimum length 100 amino acids)
                orfs = self.find_orfs(seq, min_length=300)
                
                for i, orf in enumerate(orfs):
                    gene_data = {
                        'gene_id': f"{gene_id}_ORF{i+1}",
                        'description': description,
                        'sequence': orf['sequence'],
                        'start': orf['start'],
                        'end': orf['end'],
                        'strand': orf['strand'],
                        'contig': gene_id
                    }
                    genes.append(gene_data)
            else:
                # Check annotations
                if annotations and gene_id in annotations:
                    ann = annotations[gene_id]
                    if self.is_uncharacterized(ann.get('description', '')):
                        # Extract sequence from annotation coordinates
                        start = ann['start'] - 1  # Convert to 0-based
                        end = ann['end']
                        if ann['strand'] == '+':
                            seq = str(record.seq[start:end])
                        else:
                            seq = str(record.seq[start:end].reverse_complement())
                        
                        gene_data = {
                            'gene_id': gene_id,
                            'description': ann['description'],
                            'sequence': seq,
                            'start': start,
                            'end': end,
                            'strand': ann['strand'],
                            'contig': record.id
                        }
                        genes.append(gene_data)
        
        return genes
    
    def find_orfs(self, sequence: str, min_length: int = 300) -> List[Dict]:
        """
        Find Open Reading Frames (ORFs) in a DNA sequence.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of ORF dictionaries
        """
        orfs = []
        seq = Seq(sequence)
        
        # Standard genetic code start codons
        start_codons = ['ATG', 'GTG', 'TTG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        # Search in forward strand
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3].upper()
                if codon in start_codons:
                    # Look for stop codon
                    for j in range(i+3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3].upper()
                        if stop_codon in stop_codons:
                            orf_length = j - i
                            if orf_length >= min_length:
                                orf_seq = sequence[i:j+3]
                                orfs.append({
                                    'sequence': orf_seq,
                                    'start': i,
                                    'end': j+3,
                                    'strand': '+',
                                    'frame': frame
                                })
                            break
        
        # Search in reverse complement
        rev_seq = str(seq.reverse_complement())
        for frame in range(3):
            for i in range(frame, len(rev_seq) - 2, 3):
                codon = rev_seq[i:i+3].upper()
                if codon in start_codons:
                    for j in range(i+3, len(rev_seq) - 2, 3):
                        stop_codon = rev_seq[j:j+3].upper()
                        if stop_codon in stop_codons:
                            orf_length = j - i
                            if orf_length >= min_length:
                                orf_seq = rev_seq[i:j+3]
                                orfs.append({
                                    'sequence': orf_seq,
                                    'start': len(sequence) - (j+3),
                                    'end': len(sequence) - i,
                                    'strand': '-',
                                    'frame': frame
                                })
                            break
        
        return orfs
    
    def translate_sequence(self, dna_sequence: str) -> str:
        """
        Translate DNA sequence to protein.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Protein sequence
        """
        seq = Seq(dna_sequence)
        try:
            protein = seq.translate(to_stop=True)
            return str(protein)
        except:
            return ""
    
    def blast_search(self, protein_sequence: str, database: str = "nr",
                    max_results: int = 5) -> List[Dict]:
        """
        Perform BLAST search against NCBI database.
        
        Args:
            protein_sequence: Protein sequence to search
            database: BLAST database (nr, swissprot, etc.)
            max_results: Maximum number of results to return
            
        Returns:
            List of BLAST hit dictionaries
        """
        print(f"  Performing BLAST search (this may take a while)...")
        hits = []
        
        try:
            # Use NCBI BLAST web service
            result_handle = NCBIWWW.qblast("blastp", database, protein_sequence,
                                          hitlist_size=max_results)
            blast_record = NCBIXML.read(result_handle)
            
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hits.append({
                        'title': alignment.title,
                        'accession': alignment.accession,
                        'evalue': hsp.expect,
                        'identity': hsp.identities / hsp.align_length * 100,
                        'query_coverage': hsp.align_length / len(protein_sequence) * 100,
                        'description': alignment.title.split('|')[-1] if '|' in alignment.title else alignment.title
                    })
                    break  # Only take best HSP per alignment
        except Exception as e:
            print(f"  Warning: BLAST search failed: {e}")
            print("  You may need to install BLAST+ locally or check internet connection")
        
        return hits
    
    def search_interpro(self, protein_sequence: str) -> List[Dict]:
        """
        Search InterPro for protein domains and families.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            List of domain/family information
        """
        domains = []
        
        try:
            # InterProScan API (simplified - would need actual API key for production)
            # For now, we'll use a mock approach
            print("  Searching InterPro for domains...")
            # In production, you would use InterProScan API or local installation
            # This is a placeholder
        except Exception as e:
            print(f"  Warning: InterPro search not available: {e}")
        
        return domains
    
    def analyze_gene(self, gene_data: Dict) -> Dict:
        """
        Analyze a single uncharacterized gene.
        
        Args:
            gene_data: Dictionary containing gene information
            
        Returns:
            Analysis results dictionary
        """
        print(f"\nAnalyzing gene: {gene_data['gene_id']}")
        print(f"  Description: {gene_data['description']}")
        
        # Translate to protein
        protein_seq = self.translate_sequence(gene_data['sequence'])
        if not protein_seq or len(protein_seq) < 30:
            print("  Warning: Protein sequence too short or invalid")
            return None
        
        print(f"  Protein length: {len(protein_seq)} amino acids")
        
        # Perform BLAST search
        blast_hits = self.blast_search(protein_seq, max_results=3)
        
        # Search for domains
        domains = self.search_interpro(protein_seq)
        
        # Compile results
        analysis = {
            'gene_id': gene_data['gene_id'],
            'description': gene_data['description'],
            'protein_length': len(protein_seq),
            'protein_sequence': protein_seq[:50] + "..." if len(protein_seq) > 50 else protein_seq,
            'blast_hits': blast_hits,
            'domains': domains,
            'predicted_function': self.predict_function(blast_hits, domains)
        }
        
        return analysis
    
    def predict_function(self, blast_hits: List[Dict], domains: List[Dict]) -> str:
        """
        Predict protein function based on BLAST hits and domains.
        
        Args:
            blast_hits: BLAST search results
            domains: Domain information
            
        Returns:
            Predicted function description
        """
        if not blast_hits:
            return "No significant homology found. Function remains unknown."
        
        # Get best hit
        best_hit = blast_hits[0]
        
        if best_hit['identity'] > 80 and best_hit['evalue'] < 1e-50:
            confidence = "High"
            function = f"Strongly similar to {best_hit['description']} ({best_hit['identity']:.1f}% identity)"
        elif best_hit['identity'] > 50 and best_hit['evalue'] < 1e-20:
            confidence = "Moderate"
            function = f"Moderately similar to {best_hit['description']} ({best_hit['identity']:.1f}% identity)"
        else:
            confidence = "Low"
            function = f"Weak similarity to {best_hit['description']} ({best_hit['identity']:.1f}% identity)"
        
        # Add domain information if available
        if domains:
            domain_names = [d.get('name', '') for d in domains]
            function += f". Contains domains: {', '.join(domain_names)}"
        
        return f"[{confidence} confidence] {function}"
    
    def generate_report(self, analyses: List[Dict], output_file: str):
        """
        Generate a comprehensive report of the analysis.
        
        Args:
            analyses: List of analysis results
            output_file: Output file path
        """
        print(f"\n{'='*80}")
        print("GENERATING REPORT")
        print(f"{'='*80}\n")
        
        report_lines = []
        report_lines.append("="*80)
        report_lines.append("UNCHARACTERIZED GENE ANALYSIS REPORT")
        report_lines.append("="*80)
        report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Total uncharacterized genes analyzed: {len(analyses)}")
        report_lines.append("="*80)
        report_lines.append("")
        
        for i, analysis in enumerate(analyses, 1):
            if not analysis:
                continue
                
            report_lines.append(f"\n{'='*80}")
            report_lines.append(f"GENE {i}: {analysis['gene_id']}")
            report_lines.append(f"{'='*80}")
            report_lines.append(f"Original Description: {analysis['description']}")
            report_lines.append(f"Protein Length: {analysis['protein_length']} amino acids")
            report_lines.append(f"Protein Sequence (first 50 aa): {analysis['protein_sequence']}")
            report_lines.append("")
            
            report_lines.append("PREDICTED FUNCTION:")
            report_lines.append(f"  {analysis['predicted_function']}")
            report_lines.append("")
            
            if analysis['blast_hits']:
                report_lines.append("BLAST HOMOLOGY RESULTS:")
                for j, hit in enumerate(analysis['blast_hits'], 1):
                    report_lines.append(f"  Hit {j}:")
                    report_lines.append(f"    Description: {hit['description']}")
                    report_lines.append(f"    Accession: {hit.get('accession', 'N/A')}")
                    report_lines.append(f"    Identity: {hit['identity']:.2f}%")
                    report_lines.append(f"    Query Coverage: {hit['query_coverage']:.2f}%")
                    report_lines.append(f"    E-value: {hit['evalue']:.2e}")
                    report_lines.append("")
            
            if analysis['domains']:
                report_lines.append("PROTEIN DOMAINS:")
                for domain in analysis['domains']:
                    report_lines.append(f"  - {domain.get('name', 'Unknown')}: {domain.get('description', '')}")
                report_lines.append("")
        
        report_text = "\n".join(report_lines)
        
        # Write to file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(report_text)
        
        # Also print summary
        print(report_text)
        print(f"\nReport saved to: {output_file}")
        
        # Save JSON version
        json_file = output_file.replace('.txt', '.json')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(analyses, f, indent=2, ensure_ascii=False)
        print(f"JSON report saved to: {json_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze uncharacterized genes in a genome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze genome with FASTA file
  python genome_uncharacterized_analyzer.py genome.fasta
  
  # With annotation file
  python genome_uncharacterized_analyzer.py genome.fasta -a annotations.gff
  
  # Specify email for NCBI
  python genome_uncharacterized_analyzer.py genome.fasta -e your.email@example.com
        """
    )
    
    parser.add_argument('genome_file', help='Path to genome file (FASTA format)')
    parser.add_argument('-a', '--annotations', help='Path to annotation file (GFF/GTF format)')
    parser.add_argument('-o', '--output', default='uncharacterized_genes_report.txt',
                       help='Output report file (default: uncharacterized_genes_report.txt)')
    parser.add_argument('-e', '--email', default='your.email@example.com',
                       help='Email for NCBI Entrez (required for BLAST searches)')
    parser.add_argument('-k', '--api-key', help='NCBI API key (optional, for higher rate limits)')
    parser.add_argument('-f', '--format', default='fasta',
                       choices=['fasta', 'genbank', 'embl'],
                       help='Genome file format (default: fasta)')
    parser.add_argument('--min-length', type=int, default=300,
                       help='Minimum ORF length in nucleotides (default: 300)')
    parser.add_argument('--max-genes', type=int, default=10,
                       help='Maximum number of genes to analyze (default: 10)')
    
    args = parser.parse_args()
    
    # Check if genome file exists
    if not os.path.exists(args.genome_file):
        print(f"ERROR: Genome file not found: {args.genome_file}")
        sys.exit(1)
    
    # Initialize analyzer
    print("Initializing Uncharacterized Gene Analyzer...")
    analyzer = UncharacterizedGeneAnalyzer(email=args.email, api_key=args.api_key)
    
    # Load genome
    records = analyzer.load_genome(args.genome_file, args.format)
    
    # Load annotations if provided
    annotations = analyzer.load_annotations(args.annotations)
    
    # Extract uncharacterized genes
    print("\nExtracting uncharacterized genes...")
    genes = analyzer.extract_genes_from_genome(records, annotations)
    
    if not genes:
        print("No uncharacterized genes found!")
        sys.exit(0)
    
    print(f"Found {len(genes)} uncharacterized genes")
    
    # Limit number of genes to analyze
    if len(genes) > args.max_genes:
        print(f"Limiting analysis to first {args.max_genes} genes")
        genes = genes[:args.max_genes]
    
    # Analyze each gene
    analyses = []
    for gene in genes:
        analysis = analyzer.analyze_gene(gene)
        if analysis:
            analyses.append(analysis)
    
    # Generate report
    if analyses:
        analyzer.generate_report(analyses, args.output)
    else:
        print("No genes could be analyzed.")


if __name__ == "__main__":
    main()

