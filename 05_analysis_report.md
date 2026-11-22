# DNA Sequence Analysis Report

- Run date: 2026-11-22T18:37:00+00:00
- Input file: `data/sample.fasta`
- Number of sequences: 1
- Total sequence length: 80 bp
- Mean GC content: 52.50%
- Minimum ORF length: 10 amino acids
- Candidate ORFs detected: 1

## Sequence statistics

| Sequence | Length (bp) | GC (%) | Ambiguous bases |
|---|---:|---:|---:|
| sample_gene_like_sequence Synthetic teaching sequence for pipeline demonstration | 80 | 52.50 | 0 |

## Interpretation

GC content describes the proportion of G and C bases among known A/C/G/T bases. ORFs are candidate coding regions beginning with ATG and ending at an in-frame stop codon. ORF detection alone does not prove that a region is a functional gene; biological validation and additional evidence are required.

## Output files

- `01_cleaned_sequences.fasta`
- `02_sequence_statistics.tsv`
- `03_orf_candidates.tsv`
- `04_predicted_proteins.fasta`
- `05_analysis_report.md`
