"""Prepare alignment for GEMME/ESCOTT.

This script prepares a FASTA-like alignment for use with GEMME/ESCOTT by:
- Removing lowercase letters from sequences.
- Converting gap-like characters (`.`) into `-`.
- Removing positions containing gaps in the query sequence.
- Shortening headers to a maximum length of 50 characters.

The processed alignment is saved in the same directory as the input file.
The output file will have the same name of the input `alignment_file`,
with `_processed.fasta` appended.
"""

import argparse
import string
from pathlib import Path


class FastaDict:
    """A dictionary-like object for FASTA records and sequences."""

    def __init__(self, records: list[str], sequences: list[str]):
        self.records = records
        self.sequences = sequences

    @classmethod
    def parse(cls, input_file: Path):
        """Parses a FASTA-like file into a structured dictionary.

        - Removes lowercase letters from sequences.
        - Converts gap-like characters (`.`) into `-`.
        """
        records, sequences = [], []
        translator = str.maketrans("", "", string.ascii_lowercase)
        with open(input_file, "r") as f:
            record, sequence = None, []
            for line in map(str.strip, f):
                if not line or line.startswith("#"):  # Skip empty/comment lines
                    continue
                if line.startswith(">"):
                    if record is not None:
                        records.append(record)
                        sequences.append("".join(sequence))
                    record, sequence = line.split()[0][1:], []
                else:
                    sequence.append(line.translate(translator).replace(".", "-"))
        if record is not None:
            records.append(record)
            sequences.append("".join(sequence))

        return cls(records, sequences)

    def write(self, output_file: Path):
        """Writes a FastaDict to a FASTA file."""
        with open(output_file, "w") as f:
            for record, sequence in zip(self.records, self.sequences):
                f.write(f">{record}\n{sequence}\n")


def _shorten_headers(fasta: FastaDict, max_length: int = 50) -> FastaDict:
    """Shortens FASTA headers to a maximum length."""
    return FastaDict(
        records=[record[:max_length] for record in fasta.records],
        sequences=fasta.sequences,
    )


def _remove_query_gaps(alignment: FastaDict) -> FastaDict:
    """Removes positions containing gaps in the query sequence."""
    query_sequence = alignment.sequences[0]
    gap_indices = {i for i, c in enumerate(query_sequence) if c == "-"}

    new_sequences = [
        "".join(c for i, c in enumerate(seq) if i not in gap_indices)
        for seq in alignment.sequences
    ]

    return FastaDict(alignment.records, new_sequences)


def prepare_alignment(alignment_file: Path) -> None:
    """Converts a FASTA-like alignment into FASTA format for GEMME/ESCOTT:
    - Removes lowercase letters from sequences.
    - Converts gap-like characters (`.`) into `-`.
    - Removes positions containing gaps in the query sequence.
    - Shortens headers to a maximum length of 50 characters.

    The processed alignment is saved in the same directory as the input file.
    The output file will have the same name of the input `alignment_file`,
    with `_processed.fasta` appended.

    Args:
        alignment_file: Path to the FASTA-like alignment file.

    Returns:
        None

    """
    alignment_file = Path(alignment_file)
    alignment = FastaDict.parse(alignment_file)

    # Ensure all sequences are of equal length
    seq_length = len(alignment.sequences[0])
    if any(len(seq) != seq_length for seq in alignment.sequences):
        raise ValueError("Sequences in the provided alignment have different lengths.")

    alignment = _remove_query_gaps(alignment)
    alignment = _shorten_headers(alignment)

    output_file = alignment_file.with_name(f"{alignment_file.stem}_processed.fasta")
    alignment.write(output_file)
    print(f"Alignment saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("alignment_file", help="Path to the FASTA-like alignment file.")
    args = parser.parse_args()

    prepare_alignment(args.alignment_file)
