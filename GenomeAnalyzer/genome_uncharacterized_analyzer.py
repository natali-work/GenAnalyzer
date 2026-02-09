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
import time

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez
    from Bio import ExPASy
    from Bio import SwissProt
except ImportError:
    print("ERROR: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

try:
    import requests
except ImportError:
    print("ERROR: requests is required. Install with: pip install requests")


def vprint(message: str, end: str = "\n"):
    """Verbose print with immediate flush for real-time output."""
    print(message, end=end, flush=True)


def timestamp():
    """Return current timestamp string."""
    return datetime.now().strftime("%H:%M:%S")


class UncharacterizedGeneAnalyzer:
    """Analyzes uncharacterized genes in a genome."""
    
    # Keywords that indicate uncharacterized/hypothetical proteins
    UNCHARACTERIZED_KEYWORDS = [
        'hypothetical', 'uncharacterized', 'unknown', 'putative',
        'predicted', 'similar to', 'conserved', 'domain-containing',
        'orf', 'open reading frame', 'unnamed', 'unnamed protein'
    ]
    
    # Keywords that indicate an uncharacterized BLAST hit (to be dismissed
    # when looking for homologous genes with known function)
    UNCHARACTERIZED_HIT_KEYWORDS = [
        'uncharacterized protein', 'hypothetical protein'
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
        vprint(f"  [{timestamp()}] Loading genome from {genome_file}...")
        try:
            records = list(SeqIO.parse(genome_file, file_format))
            vprint(f"  [{timestamp()}] Loaded {len(records)} sequences")
            return records
        except Exception as e:
            vprint(f"  [{timestamp()}] Error loading genome: {e}")
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
            vprint(f"  [{timestamp()}] No annotation file provided. Will extract from sequence headers.")
            return annotations
            
        vprint(f"  [{timestamp()}] Loading annotations from {annotation_file}...")
        try:
            with open(annotation_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9:
                        # Only process CDS entries (coding sequences)
                        if parts[2] != 'CDS':
                            continue
                        attributes = parts[8]
                        # Extract gene ID and description
                        gene_id = None
                        description = ""
                        
                        # Parse GTF format: key "value"; key "value";
                        # or GFF format: key=value;key=value
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if not attr:
                                continue
                            
                            # GTF format: gene_id "value"
                            if ' "' in attr:
                                key_value = attr.split(' "', 1)
                                if len(key_value) == 2:
                                    key = key_value[0].strip()
                                    value = key_value[1].rstrip('"').strip()
                                    if key in ('gene_id', 'locus_tag'):
                                        gene_id = value
                                    if key == 'product':
                                        description = value
                            # GFF format: key=value
                            elif '=' in attr:
                                key_value = attr.split('=', 1)
                                if len(key_value) == 2:
                                    key = key_value[0].strip()
                                    value = key_value[1].strip()
                                    if key in ('ID', 'gene_id', 'locus_tag'):
                                        gene_id = value
                                    if key in ('Name', 'gene', 'product'):
                                        description = value
                        
                        if gene_id and gene_id not in annotations:
                            annotations[gene_id] = {
                                'description': description,
                                'type': parts[2],
                                'start': int(parts[3]),
                                'end': int(parts[4]),
                                'strand': parts[6],
                                'contig': parts[0]
                            }
            vprint(f"  [{timestamp()}] Loaded {len(annotations)} annotations")
        except Exception as e:
            vprint(f"  [{timestamp()}] Warning: Could not load annotations: {e}")
            
        return annotations
    
    def is_characterized_hit(self, hit_description: str) -> bool:
        """
        Check if a BLAST hit describes a characterized (known-function) protein.
        
        Hits matching 'uncharacterized protein' or 'hypothetical protein' are
        dismissed so the algorithm can keep scanning for homology to a gene
        with a known function.
        
        Args:
            hit_description: BLAST hit description/title
            
        Returns:
            True if the hit is characterized (i.e. NOT hypothetical/uncharacterized)
        """
        desc_lower = hit_description.lower()
        return not any(kw in desc_lower for kw in self.UNCHARACTERIZED_HIT_KEYWORDS)
    
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
        
        # Build a lookup of records by ID
        record_lookup = {record.id: record for record in records}
        
        # If we have annotations, use them to find uncharacterized genes
        if annotations:
            for gene_id, ann in annotations.items():
                # Check if this gene is uncharacterized
                if self.is_uncharacterized(ann.get('description', '')):
                    # Find the corresponding sequence record
                    contig_id = ann.get('contig', '')
                    record = record_lookup.get(contig_id)
                    
                    if record:
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
                            'contig': contig_id
                        }
                        genes.append(gene_data)
        else:
            # No annotations - check sequence headers
            for record in records:
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
        
        Hits whose descriptions match uncharacterized/hypothetical protein are
        automatically filtered out so that only homology to characterized
        (known-function) proteins is reported.  To compensate for the
        filtering, extra hits are requested from NCBI.
        
        Args:
            protein_sequence: Protein sequence to search
            database: BLAST database (nr, swissprot, etc.)
            max_results: Maximum number of characterized results to return
            
        Returns:
            List of BLAST hit dictionaries (only characterized hits)
        """
        # Request extra hits to compensate for filtering out uncharacterized ones
        request_size = max_results * 5
        
        vprint(f"  [{timestamp()}] Starting BLAST search against '{database}' database...")
        vprint(f"  [{timestamp()}] Sequence length: {len(protein_sequence)} amino acids")
        vprint(f"  [{timestamp()}] Requesting up to {request_size} hits (will keep top {max_results} characterized)")
        vprint(f"  [{timestamp()}] Submitting query to NCBI servers (this typically takes 1-5 minutes)...")
        hits = []
        
        start_time = time.time()
        
        try:
            # Use NCBI BLAST web service
            vprint(f"  [{timestamp()}] Waiting for NCBI BLAST response", end="")
            result_handle = NCBIWWW.qblast("blastp", database, protein_sequence,
                                          hitlist_size=request_size)
            
            elapsed = time.time() - start_time
            vprint(f"\n  [{timestamp()}] BLAST query completed in {elapsed:.1f} seconds")
            vprint(f"  [{timestamp()}] Parsing BLAST results...")
            
            blast_record = NCBIXML.read(result_handle)
            
            vprint(f"  [{timestamp()}] Found {len(blast_record.alignments)} alignments")
            
            skipped = 0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    description = alignment.title.split('|')[-1] if '|' in alignment.title else alignment.title
                    
                    # Skip uncharacterized / hypothetical protein hits
                    if not self.is_characterized_hit(description):
                        skipped += 1
                        break  # Skip to next alignment
                    
                    hits.append({
                        'title': alignment.title,
                        'accession': alignment.accession,
                        'evalue': hsp.expect,
                        'identity': hsp.identities / hsp.align_length * 100,
                        'query_coverage': hsp.align_length / len(protein_sequence) * 100,
                        'description': description
                    })
                    break  # Only take best HSP per alignment
                
                # Stop once we have enough characterized hits
                if len(hits) >= max_results:
                    break
            
            vprint(f"  [{timestamp()}] Processed {len(hits)} characterized BLAST hits (skipped {skipped} uncharacterized/hypothetical)")
            
        except Exception as e:
            elapsed = time.time() - start_time
            vprint(f"\n  [{timestamp()}] WARNING: BLAST search failed after {elapsed:.1f}s: {e}")
            vprint(f"  [{timestamp()}] You may need to install BLAST+ locally or check internet connection")
        
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
            vprint(f"  [{timestamp()}] InterPro domain search (placeholder - not implemented)")
            # In production, you would use InterProScan API or local installation
            # This is a placeholder
        except Exception as e:
            vprint(f"  [{timestamp()}] Warning: InterPro search not available: {e}")
        
        return domains
    
    def analyze_gene(self, gene_data: Dict) -> Dict:
        """
        Analyze a single uncharacterized gene.
        
        Args:
            gene_data: Dictionary containing gene information
            
        Returns:
            Analysis results dictionary
        """
        gene_start_time = time.time()
        
        vprint(f"\n{'-'*60}")
        vprint(f"[{timestamp()}] ANALYZING GENE: {gene_data['gene_id']}")
        vprint(f"{'-'*60}")
        vprint(f"  [{timestamp()}] Description: {gene_data['description'][:80]}...")
        vprint(f"  [{timestamp()}] Location: {gene_data.get('contig', 'N/A')} : {gene_data.get('start', 'N/A')}-{gene_data.get('end', 'N/A')} ({gene_data.get('strand', 'N/A')})")
        vprint(f"  [{timestamp()}] Nucleotide sequence length: {len(gene_data['sequence'])} bp")
        
        # Translate to protein
        vprint(f"  [{timestamp()}] Translating DNA to protein sequence...")
        protein_seq = self.translate_sequence(gene_data['sequence'])
        if not protein_seq or len(protein_seq) < 30:
            vprint(f"  [{timestamp()}] WARNING: Protein sequence too short or invalid (length: {len(protein_seq) if protein_seq else 0})")
            vprint(f"  [{timestamp()}] Skipping this gene.")
            return None
        
        vprint(f"  [{timestamp()}] Protein length: {len(protein_seq)} amino acids")
        
        # Perform BLAST search
        vprint(f"  [{timestamp()}] Starting homology search...")
        blast_hits = self.blast_search(protein_seq, max_results=3)
        
        # Search for domains
        vprint(f"  [{timestamp()}] Searching for protein domains...")
        domains = self.search_interpro(protein_seq)
        
        # Compile results
        vprint(f"  [{timestamp()}] Compiling analysis results...")
        analysis = {
            'gene_id': gene_data['gene_id'],
            'description': gene_data['description'],
            'protein_length': len(protein_seq),
            'protein_sequence': protein_seq[:50] + "..." if len(protein_seq) > 50 else protein_seq,
            'blast_hits': blast_hits,
            'domains': domains,
            'predicted_function': self.predict_function(blast_hits, domains)
        }
        
        gene_elapsed = time.time() - gene_start_time
        vprint(f"  [{timestamp()}] Gene analysis completed in {gene_elapsed:.1f} seconds")
        vprint(f"  [{timestamp()}] Predicted function: {analysis['predicted_function'][:70]}...")
        
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
        vprint(f"\n{'='*80}")
        vprint(f"[{timestamp()}] GENERATING REPORT")
        vprint(f"{'='*80}\n")
        
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
        vprint(report_text)
        vprint(f"\n[{timestamp()}] Report saved to: {output_file}")
        
        # Save JSON version
        json_file = output_file.replace('.txt', '.json')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(analyses, f, indent=2, ensure_ascii=False)
        vprint(f"[{timestamp()}] JSON report saved to: {json_file}")
    
    def generate_summary_table(self, genes: List[Dict], output_file: str):
        """
        Generate a summary table of all uncharacterized genes (no BLAST).
        
        Args:
            genes: List of gene dictionaries
            output_file: Output CSV file path
        """
        vprint(f"\n{'='*80}")
        vprint(f"[{timestamp()}] GENERATING SUMMARY TABLE")
        vprint(f"{'='*80}\n")
        
        # CSV output
        csv_lines = []
        csv_lines.append("Gene_ID,Description,Contig,Start,End,Strand,Nucleotide_Length,Protein_Length")
        
        for gene in genes:
            protein_seq = self.translate_sequence(gene['sequence'])
            protein_len = len(protein_seq) if protein_seq else 0
            nuc_len = len(gene['sequence'])
            
            # Escape commas in description
            description = gene['description'].replace(',', ';')
            
            csv_lines.append(f"{gene['gene_id']},{description},{gene['contig']},{gene['start']},{gene['end']},{gene['strand']},{nuc_len},{protein_len}")
        
        csv_text = "\n".join(csv_lines)
        
        # Write CSV
        csv_file = output_file.replace('.txt', '_summary.csv')
        with open(csv_file, 'w', encoding='utf-8') as f:
            f.write(csv_text)
        vprint(f"[{timestamp()}] Summary CSV saved to: {csv_file}")
        
        # Also generate a formatted text table
        table_lines = []
        table_lines.append("="*120)
        table_lines.append("UNCHARACTERIZED GENES SUMMARY TABLE")
        table_lines.append("="*120)
        table_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        table_lines.append(f"Total uncharacterized genes: {len(genes)}")
        table_lines.append("="*120)
        table_lines.append("")
        table_lines.append(f"{'No.':<6} {'Gene ID':<20} {'Protein (aa)':<12} {'Strand':<8} {'Start':<12} {'End':<12} {'Description':<50}")
        table_lines.append("-"*120)
        
        for i, gene in enumerate(genes, 1):
            protein_seq = self.translate_sequence(gene['sequence'])
            protein_len = len(protein_seq) if protein_seq else 0
            desc = gene['description'][:47] + "..." if len(gene['description']) > 50 else gene['description']
            table_lines.append(f"{i:<6} {gene['gene_id']:<20} {protein_len:<12} {gene['strand']:<8} {gene['start']:<12} {gene['end']:<12} {desc:<50}")
        
        table_lines.append("-"*120)
        table_lines.append(f"Total: {len(genes)} uncharacterized genes")
        
        table_text = "\n".join(table_lines)
        
        # Write text table
        table_file = output_file.replace('.txt', '_summary.txt')
        with open(table_file, 'w', encoding='utf-8') as f:
            f.write(table_text)
        
        vprint(table_text)
        vprint(f"\n[{timestamp()}] Summary table saved to: {table_file}")
        vprint(f"[{timestamp()}] Summary CSV saved to: {csv_file}")
        
        return csv_file, table_file


def main():
    parser = argparse.ArgumentParser(
        description="Analyze uncharacterized genes in a genome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze genome with FASTA file (interactive batch mode, 20 genes at a time)
  python genome_uncharacterized_analyzer.py genome.fasta
  
  # With annotation file
  python genome_uncharacterized_analyzer.py genome.fasta -a annotations.gff
  
  # Specify email for NCBI
  python genome_uncharacterized_analyzer.py genome.fasta -e your.email@example.com
  
  # Analyze in batches of 50 genes
  python genome_uncharacterized_analyzer.py genome.fasta --batch-size 50
  
  # Generate summary only (no BLAST, fast)
  python genome_uncharacterized_analyzer.py genome.fasta --summary-only
  
  # Auto-continue through all batches (results saved after each batch)
  python genome_uncharacterized_analyzer.py genome.fasta --auto-continue
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
    parser.add_argument('--summary-only', action='store_true',
                       help='Generate summary table of ALL genes without BLAST (fast)')
    parser.add_argument('--batch-size', type=int, default=20,
                       help='Number of genes to analyze per batch (default: 20)')
    parser.add_argument('--auto-continue', action='store_true',
                       help='Automatically continue to next batch without prompting (results saved after each batch)')
    
    args = parser.parse_args()
    
    # Check if genome file exists
    if not os.path.exists(args.genome_file):
        vprint(f"[{timestamp()}] ERROR: Genome file not found: {args.genome_file}")
        sys.exit(1)
    
    # Initialize analyzer
    vprint(f"[{timestamp()}] Initializing Uncharacterized Gene Analyzer...")
    analyzer = UncharacterizedGeneAnalyzer(email=args.email, api_key=args.api_key)
    
    # Load genome
    vprint(f"[{timestamp()}] Loading genome file...")
    records = analyzer.load_genome(args.genome_file, args.format)
    
    # Load annotations if provided
    vprint(f"[{timestamp()}] Loading annotations...")
    annotations = analyzer.load_annotations(args.annotations)
    
    # Extract uncharacterized genes
    vprint(f"\n[{timestamp()}] Extracting uncharacterized genes...")
    genes = analyzer.extract_genes_from_genome(records, annotations)
    
    if not genes:
        vprint(f"[{timestamp()}] No uncharacterized genes found!")
        sys.exit(0)
    
    vprint(f"[{timestamp()}] Found {len(genes)} uncharacterized genes")
    
    # Always generate summary table first
    vprint(f"\n[{timestamp()}] Generating summary table of ALL uncharacterized genes...")
    analyzer.generate_summary_table(genes, args.output)
    
    # If summary-only mode, stop here
    if args.summary_only:
        vprint(f"\n[{timestamp()}] [Summary-only mode] Skipping BLAST analysis.")
        vprint(f"[{timestamp()}] To run BLAST analysis, remove --summary-only flag.")
        sys.exit(0)
    
    # Batch processing with user interaction
    batch_size = args.batch_size
    total_genes = len(genes)
    all_analyses = []
    current_index = 0
    batch_number = 0
    total_start_time = time.time()
    
    vprint(f"\n{'='*80}")
    vprint(f"[{timestamp()}] STARTING BLAST ANALYSIS")
    vprint(f"{'='*80}")
    vprint(f"[{timestamp()}] Total uncharacterized genes to analyze: {total_genes}")
    vprint(f"[{timestamp()}] Batch size: {batch_size} genes per batch")
    vprint(f"[{timestamp()}] Total batches: {(total_genes + batch_size - 1) // batch_size}")
    vprint(f"[{timestamp()}] NOTE: Each BLAST search takes 1-5 minutes. A batch of {batch_size} genes may take 20-100 minutes.")
    vprint(f"{'='*80}")
    
    while current_index < total_genes:
        batch_number += 1
        batch_end = min(current_index + batch_size, total_genes)
        batch_genes = genes[current_index:batch_end]
        batch_start_time = time.time()
        
        vprint(f"\n{'='*80}")
        vprint(f"[{timestamp()}] BATCH {batch_number}: Analyzing genes {current_index + 1} to {batch_end} of {total_genes}")
        vprint(f"{'='*80}")
        
        # Analyze each gene in this batch with BLAST
        batch_analyses = []
        for i, gene in enumerate(batch_genes, 1):
            overall_index = current_index + i
            vprint(f"\n[{timestamp()}] ============================================================")
            vprint(f"[{timestamp()}] GENE {overall_index}/{total_genes} (Batch {batch_number}, Item {i}/{len(batch_genes)})")
            vprint(f"[{timestamp()}] ============================================================")
            analysis = analyzer.analyze_gene(gene)
            if analysis:
                batch_analyses.append(analysis)
        
        all_analyses.extend(batch_analyses)
        batch_elapsed = time.time() - batch_start_time
        total_elapsed = time.time() - total_start_time
        
        # Display batch summary
        vprint(f"\n{'='*80}")
        vprint(f"[{timestamp()}] BATCH {batch_number} SUMMARY")
        vprint(f"{'='*80}")
        vprint(f"[{timestamp()}] Genes analyzed in this batch: {len(batch_genes)}")
        vprint(f"[{timestamp()}] Successful analyses: {len(batch_analyses)}")
        vprint(f"[{timestamp()}] Batch time: {batch_elapsed/60:.1f} minutes ({batch_elapsed:.0f} seconds)")
        vprint(f"[{timestamp()}] Total time so far: {total_elapsed/60:.1f} minutes")
        vprint(f"[{timestamp()}] Total genes analyzed so far: {batch_end}")
        vprint(f"[{timestamp()}] Remaining genes: {total_genes - batch_end}")
        vprint("")
        
        if batch_analyses:
            vprint("Genes in this batch:")
            vprint("-" * 80)
            for analysis in batch_analyses:
                func_preview = analysis['predicted_function'][:60] + "..." if len(analysis['predicted_function']) > 60 else analysis['predicted_function']
                vprint(f"  * {analysis['gene_id']}: {func_preview}")
            vprint("-" * 80)
        
        # Save results after each batch (cumulative)
        if all_analyses:
            vprint(f"\n[{timestamp()}] Saving cumulative results (all batches so far)...")
            analyzer.generate_report(all_analyses, args.output)
            vprint(f"[{timestamp()}] Results saved to: {args.output}")
        
        # Update current index
        current_index = batch_end
        
        # Check if there are more genes to process
        if current_index < total_genes:
            remaining = total_genes - current_index
            remaining_batches = (remaining + batch_size - 1) // batch_size
            avg_time_per_gene = total_elapsed / batch_end if batch_end > 0 else 0
            estimated_remaining = avg_time_per_gene * remaining
            
            vprint(f"\n[{timestamp()}] {remaining} genes remaining ({remaining_batches} more batch(es))")
            vprint(f"[{timestamp()}] Estimated time for remaining genes: {estimated_remaining/60:.1f} minutes")
            vprint("")
            
            # Auto-continue or prompt user
            if args.auto_continue:
                vprint(f"[{timestamp()}] [Auto-continue mode] Automatically proceeding to next batch...")
                vprint(f"[{timestamp()}] Results have been saved. Continuing...")
            else:
                # Prompt user for next action
                user_input = None
                while True:
                    try:
                        user_input = input(f"[{timestamp()}] Would you like to analyze the next batch? (y/n/save): ").strip().lower()
                    except (EOFError, KeyboardInterrupt):
                        vprint(f"\n[{timestamp()}] Input interrupted. Stopping analysis.")
                        vprint(f"[{timestamp()}] Analyzed {current_index} of {total_genes} genes.")
                        vprint(f"[{timestamp()}] Results have been saved to: {args.output}")
                        user_input = 'n'  # Set to 'n' to trigger break
                        break
                    
                    if user_input in ('y', 'yes'):
                        vprint(f"\n[{timestamp()}] Continuing to next batch...")
                        break
                    elif user_input in ('n', 'no'):
                        vprint(f"\n[{timestamp()}] Stopping analysis.")
                        vprint(f"[{timestamp()}] Analyzed {current_index} of {total_genes} genes.")
                        break
                    elif user_input == 'save':
                        # Save current progress and continue
                        if all_analyses:
                            analyzer.generate_report(all_analyses, args.output)
                            vprint(f"\n[{timestamp()}] Results saved. Continuing...")
                        else:
                            vprint(f"[{timestamp()}] No analyses to save yet.")
                        continue
                    else:
                        vprint(f"[{timestamp()}] Please enter 'y' (yes), 'n' (no), or 'save' to save current progress.")
                
                if user_input in ('n', 'no'):
                    break
        else:
            total_elapsed = time.time() - total_start_time
            vprint(f"\n[{timestamp()}] All {total_genes} genes have been analyzed!")
            vprint(f"[{timestamp()}] Total analysis time: {total_elapsed/60:.1f} minutes")
    
    # Generate final report (final save - results already saved after each batch)
    if all_analyses:
        total_elapsed = time.time() - total_start_time
        vprint(f"\n{'='*80}")
        vprint(f"[{timestamp()}] FINAL SUMMARY")
        vprint(f"{'='*80}")
        vprint(f"[{timestamp()}] Total genes analyzed: {len(all_analyses)}")
        vprint(f"[{timestamp()}] Total time: {total_elapsed/60:.1f} minutes ({total_elapsed/3600:.2f} hours)")
        vprint(f"[{timestamp()}] Final report saved to: {args.output}")
        # Final save (already saved after each batch, but ensure it's up to date)
        analyzer.generate_report(all_analyses, args.output)
    else:
        vprint(f"[{timestamp()}] No genes could be analyzed.")


if __name__ == "__main__":
    main()

