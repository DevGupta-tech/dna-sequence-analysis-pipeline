#!/usr/bin/env python3
"""
DNA Sequence Analysis Pipeline
Runs sequence cleaning, basic statistics, six-frame ORF detection,
protein translation, and report generation in one command.
"""

from pathlib import Path
import argparse
import csv
from datetime import datetime

# --------------------------- PARAMETERS ---------------------------
INPUT_FASTA = "data/sample.fasta"
OUTPUT_DIR = "outputs"
MIN_ORF_AA_LENGTH = 10
START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
VALID_BASES = set("ACGTN")
# ------------------------------------------------------------------

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def read_fasta(path):
    """Read a FASTA file and return [(identifier, sequence), ...]."""
    records = []
    identifier = None
    chunks = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(chunks).upper()))
                identifier = line[1:].strip() or f"sequence_{len(records)+1}"
                chunks = []
            else:
                if identifier is None:
                    raise ValueError(
                        f"FASTA sequence found before a header at line {line_number}."
                    )
                chunks.append("".join(line.split()))
    if identifier is not None:
        records.append((identifier, "".join(chunks).upper()))
    if not records:
        raise ValueError("No FASTA records were found in the input file.")
    return records

def clean_and_validate(records):
    cleaned = []
    for identifier, sequence in records:
        sequence = sequence.upper().replace("U", "T")
        invalid = sorted(set(sequence) - VALID_BASES)
        if invalid:
            raise ValueError(
                f"Sequence '{identifier}' contains unsupported symbols: {invalid}"
            )
        if not sequence:
            raise ValueError(f"Sequence '{identifier}' is empty.")
        cleaned.append((identifier, sequence))
    return cleaned

def gc_percentage(sequence):
    """Calculate GC percentage, excluding ambiguous N bases."""
    known_bases = [base for base in sequence if base in "ACGT"]
    if not known_bases:
        return 0.0
    return 100.0 * sum(base in "GC" for base in known_bases) / len(known_bases)

def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def translate(sequence):
    amino_acids = []
    for index in range(0, len(sequence) - 2, 3):
        amino_acids.append(CODON_TABLE.get(sequence[index:index+3], "X"))
    return "".join(amino_acids)

def find_orfs(identifier, sequence, min_aa_length):
    """Find ATG-to-stop ORFs in all three forward and three reverse frames.

    Coordinates are 1-based and refer to the strand being scanned.
    The stop codon is excluded from the reported protein sequence.
    """
    found = []
    strands = [("+", sequence), ("-", reverse_complement(sequence))]
    for strand, oriented_sequence in strands:
        for frame in range(3):
            for start in range(frame, len(oriented_sequence) - 2, 3):
                if oriented_sequence[start:start+3] not in START_CODONS:
                    continue
                for stop in range(start + 3, len(oriented_sequence) - 2, 3):
                    if oriented_sequence[stop:stop+3] in STOP_CODONS:
                        coding_dna = oriented_sequence[start:stop]
                        protein = translate(coding_dna)
                        if len(protein) >= min_aa_length:
                            found.append({
                                "sequence_id": identifier,
                                "strand": strand,
                                "frame": frame + 1,
                                "start": start + 1,
                                "end": stop + 3,
                                "length_nt": len(coding_dna),
                                "length_aa": len(protein),
                                "dna": coding_dna,
                                "protein": protein,
                            })
                        break
    return found

def write_fasta(records, path, width=60):
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for index in range(0, len(sequence), width):
                handle.write(sequence[index:index+width] + "\n")

def run_pipeline(input_path, output_dir, min_aa_length):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Read, clean, validate, and preserve the intermediate FASTA.
    raw_records = read_fasta(input_path)
    records = clean_and_validate(raw_records)
    cleaned_path = output_dir / "01_cleaned_sequences.fasta"
    write_fasta(records, cleaned_path)

    # Step 2: Calculate sequence-level statistics.
    stats_path = output_dir / "02_sequence_statistics.tsv"
    stats = []
    for identifier, sequence in records:
        row = {
            "sequence_id": identifier,
            "length_bp": len(sequence),
            "gc_percent": round(gc_percentage(sequence), 2),
            "ambiguous_bases": sequence.count("N"),
        }
        stats.append(row)
    with stats_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(stats)

    # Step 3: Identify candidate ORFs in six reading frames.
    all_orfs = []
    for identifier, sequence in records:
        all_orfs.extend(find_orfs(identifier, sequence, min_aa_length))
    orf_path = output_dir / "03_orf_candidates.tsv"
    orf_fields = [
        "sequence_id", "strand", "frame", "start", "end",
        "length_nt", "length_aa", "dna"
    ]
    with orf_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=orf_fields, delimiter="\t")
        writer.writeheader()
        for orf in all_orfs:
            writer.writerow({key: orf[key] for key in orf_fields})

    # Step 4: Translate candidate ORFs and preserve protein FASTA.
    protein_path = output_dir / "04_predicted_proteins.fasta"
    protein_records = []
    for index, orf in enumerate(all_orfs, start=1):
        name = (
            f"ORF_{index}|{orf['sequence_id']}|strand={orf['strand']}"
            f"|frame={orf['frame']}|aa={orf['length_aa']}"
        )
        protein_records.append((name, orf["protein"]))
    write_fasta(protein_records, protein_path)

    # Step 5: Create a human-readable summary report.
    report_path = output_dir / "05_analysis_report.md"
    total_bases = sum(row["length_bp"] for row in stats)
    average_gc = sum(row["gc_percent"] for row in stats) / len(stats)
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# DNA Sequence Analysis Report\n\n")
        handle.write(f"- Run date: {datetime.now().astimezone().isoformat(timespec='seconds')}\n")
        handle.write(f"- Input file: `{input_path}`\n")
        handle.write(f"- Number of sequences: {len(records)}\n")
        handle.write(f"- Total sequence length: {total_bases} bp\n")
        handle.write(f"- Mean GC content: {average_gc:.2f}%\n")
        handle.write(f"- Minimum ORF length: {min_aa_length} amino acids\n")
        handle.write(f"- Candidate ORFs detected: {len(all_orfs)}\n\n")
        handle.write("## Sequence statistics\n\n")
        handle.write("| Sequence | Length (bp) | GC (%) | Ambiguous bases |\n")
        handle.write("|---|---:|---:|---:|\n")
        for row in stats:
            handle.write(
                f"| {row['sequence_id']} | {row['length_bp']} | "
                f"{row['gc_percent']:.2f} | {row['ambiguous_bases']} |\n"
            )
        handle.write("\n## Interpretation\n\n")
        handle.write(
            "GC content describes the proportion of G and C bases among known "
            "A/C/G/T bases. ORFs are candidate coding regions beginning with "
            "ATG and ending at an in-frame stop codon. ORF detection alone does "
            "not prove that a region is a functional gene; biological validation "
            "and additional evidence are required.\n"
        )
        handle.write("\n## Output files\n\n")
        for filename in [
            cleaned_path.name, stats_path.name, orf_path.name,
            protein_path.name, report_path.name
        ]:
            handle.write(f"- `{filename}`\n")

    print("Pipeline completed successfully.")
    print(f"Sequences analysed: {len(records)}")
    print(f"Candidate ORFs found: {len(all_orfs)}")
    print(f"Results saved in: {output_dir.resolve()}")

def main():
    parser = argparse.ArgumentParser(
        description="Run a reproducible DNA sequence analysis pipeline."
    )
    parser.add_argument("--input", default=INPUT_FASTA, help="Input FASTA file")
    parser.add_argument("--output", default=OUTPUT_DIR, help="Output directory")
    parser.add_argument(
        "--min-orf-aa", type=int, default=MIN_ORF_AA_LENGTH,
        help="Minimum ORF length in amino acids"
    )
    args = parser.parse_args()
    if args.min_orf_aa < 1:
        parser.error("--min-orf-aa must be at least 1")
    run_pipeline(Path(args.input), Path(args.output), args.min_orf_aa)

if __name__ == "__main__":
    main()
