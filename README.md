# DNA Sequence Analysis Automation Pipeline

## Project overview
This project is a small, reproducible bioinformatics workflow. It accepts DNA
sequences in FASTA format and performs several linked analysis steps using one
Python command.

## What the pipeline does
1. Reads FASTA records, converts U to T, and validates sequence symbols.
2. Saves cleaned sequences as an intermediate FASTA file.
3. Calculates sequence length, GC percentage, and ambiguous-base count.
4. Searches all six reading frames for candidate ATG-to-stop ORFs.
5. Translates candidate ORFs into amino-acid sequences.
6. Writes a Markdown report and tabular results.

**Important:** ORF prediction produces candidate coding regions; it does not
prove that a sequence is a real or functional gene. The included sequence is
synthetic teaching data, not a biological sample.

## Requirements
- Python 3.9 or newer
- No third-party Python packages are required.

## Run from a clean checkout
Open a terminal in the project folder and run:

```bash
python pipeline.py
```

The script uses `data/sample.fasta` by default and creates the `outputs/`
folder automatically.

To analyse another FASTA file:

```bash
python pipeline.py --input data/my_sequences.fasta --output results
```

To change the minimum ORF size:

```bash
python pipeline.py --min-orf-aa 20
```

## Output files
| File | Description |
|---|---|
| `01_cleaned_sequences.fasta` | Validated and normalised DNA sequences |
| `02_sequence_statistics.tsv` | Length, GC percentage, ambiguous-base count |
| `03_orf_candidates.tsv` | Candidate ORF coordinates and DNA sequences |
| `04_predicted_proteins.fasta` | Translated candidate protein sequences |
| `05_analysis_report.md` | Summary and interpretation |

Intermediate files are deliberately retained so each stage can be inspected.

## Parameters
The default input path, output folder, minimum ORF length, start codons, and
stop codons are defined together near the top of `pipeline.py`. Command-line
options can override the input, output folder, and minimum ORF length.

## Reproducibility
The project uses only the Python standard library. No external database,
network connection, API key, or manually installed package is needed. Starting
with the same input FASTA and parameters produces the same biological results.
The report also records the run date.

## Limitations
- The ORF search uses ATG as the start codon and the three standard stop codons.
- It scans six reading frames but does not model introns, promoters, gene
  regulation, or organism-specific genetic codes.
- ORFs are candidates, not confirmed genes.
- Coordinates are 1-based on the strand being scanned.

## Suggested repository structure
```text
dna-sequence-analysis-pipeline/
├── data/
│   └── sample.fasta
├── outputs/
│   └── (generated when the pipeline runs)
├── pipeline.py
├── requirements.txt
├── .gitignore
└── README.md
```
