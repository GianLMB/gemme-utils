"""Reformat PDB file such that it can be used as input for ESCOTT,
by removing all residues that are not part of the input sequence.

Requirements:
- BioPython
"""

try:
    from Bio import PDB
except ImportError:
    raise ImportError("BioPython is required to run this script.")

import argparse
from pathlib import Path


class ResidueSelect(PDB.Select):
    """Selects residues based on their ids."""

    def __init__(self, residue_numbers: list[int]):
        self.residue_numbers = residue_numbers

    def accept_residue(self, residue) -> bool:
        return residue.id[1] in self.residue_numbers


def _extract_query(file_path: str) -> str:
    """Extract lines between the first and second lines starting with '>'.

    Args:
        file_path (str): Path to the input file.

    Returns:
        list[str]: Lines between the first and second lines starting with '>'.
    """
    lines_between = []
    with open(file_path, "r") as file:
        header = False
        for line in file:
            if line.startswith(">"):
                if header:
                    break
                header = True
            elif header:
                lines_between.append(line.strip())
    return "".join(lines_between)


def prepare_pdb(pdb_file: str, alignment_file: str) -> None:
    """Extracts sequence from PDB file and save a new PDB file with only
    the residues of the input sequence.

    Args:
        pdb_file (str): Path to the input PDB file.
        alignment_file (str): Path to the alignment file in FASTA format,
            with the query sequence as first entry.

    Returns:
        None
    """
    # Get query sequence from the alignment file
    sequence = _extract_query(alignment_file)

    # Load the PDB file
    pdb_file = Path(pdb_file)
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    # Extract the sequence from the PDB file
    ppb = PDB.PPBuilder()
    for model in structure:
        for chain in model:
            chain_sequence = "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])
            if sequence in chain_sequence:
                residue_numbers = [
                    residue.id[1] for residue in chain if sequence in chain_sequence
                ]
                break
    io = PDB.PDBIO()
    io.set_structure(structure)

    output_pdb_file = pdb_file.with_name(f"{pdb_file.stem}_processed.pdb").as_posix()
    io.save(output_pdb_file, ResidueSelect(residue_numbers))
    print(f"Processed PDB file saved to {output_pdb_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pdb_file", help="Path to the input PDB file.")
    parser.add_argument(
        "alignment_file",
        help="""Path to the alignment file in FASTA format,
with the query sequence as first entry.
""",
    )
    args = parser.parse_args()

    prepare_pdb(args.pdb_file, args.alignment_file)
