"""Reformat PDB file such that it can be used as input for ESCOTT,
by removing all residues that are not part of the query sequence.

Query sequence is extracted from an alignment file in FASTA format.
The common subsequence between the query and PDB sequence is extracted,
and the PDB file is trimmed to only include these residues (output
file is saved with `_processed.pdb` appended).

If the query sequence (and thus the alignment) is longer than the PDB
sequence, the script will produce an additional alignment file
(with `_trimmed.fasta` appended) with only the common subsequence.

If an exact local match between the query and PDB sequence is not found,
the script will attempt to find the best local alignment between the two
sequences. If the alignment score is above an acceptance threshold, the
script will issue a warning and proceed with the best match, otherwise
it will raise an error.

Requirements:
- BioPython
"""

try:
    from Bio import PDB, SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.PDB.Structure import Structure
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError("BioPython is required to run this script.")

import argparse
import warnings
from pathlib import Path


# Minimum score for local alignment to consider match between query and PDB
_SCORE_THRESHOLD = 0.8


class ResidueSelect(PDB.Select):
    """Selects residues based on their ids."""

    def __init__(self, residue_numbers: list[int], chain_id: str) -> None:
        self.residue_numbers = residue_numbers
        self.chain_id = chain_id
        
    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue) -> bool:
        return residue.id[1] in self.residue_numbers


def _extract_alignment(alignment_file: Path) -> list[SeqRecord]:
    """Extracts sequences from an MSA file in FASTA format."""
    return list(SeqIO.parse(alignment_file, "fasta"))


def _trim_alignment(alignment: list[SeqRecord], trimmed_query_sequence: str) -> list[SeqRecord]:
    """Trim the alignment to match the PDB sequence."""

    if not "-" in trimmed_query_sequence:
        query_sequence = str(alignment[0].seq)
        start_idx = query_sequence.find(trimmed_query_sequence)
        end_idx = start_idx + len(trimmed_query_sequence)
        trimmed_alignment = [
            SeqRecord(seq=record.seq[start_idx:end_idx], id=record.id, description="")
            for record in alignment
        ]
    else:
        gap_indices = {i for i, c in enumerate(trimmed_query_sequence) if c == "-"}
        trimmed_alignment = [
            SeqRecord(
                seq=Seq("".join([c for i, c in enumerate(record.seq) if i not in gap_indices])),
                id=record.id,
                description="",
            )
            for record in alignment
        ]

    return trimmed_alignment


def _process_best_match(
    best_match: dict[str, any], query_sequence: str, aligner: PairwiseAligner
) -> tuple[list[int], str]:
    """Processes the best match between query and PDB sequence."""
    # get best matching chain
    chain = best_match["chain"]
    chain_sequence = best_match["chain_sequence"]
    chain = list(chain)

    # get alignment positions and gap mask
    alignment = aligner.align(query_sequence, chain_sequence)[0]
    indices = alignment.indices
    gap_mask = [(i != -1 and j != -1) for i, j in zip(indices[0], indices[1])]

    # get residue numbers
    residue_numbers = [chain[c].id[1] for c, match in zip(indices[1], gap_mask) if match]

    # get trimmed query sequence, with gaps in non-matching positions
    trimmed_query_sequence = [
        query_sequence[i] if i in indices[0] else "-" for i in range(len(query_sequence))
    ]
    # add gaps for missing pdb residues
    missing_indices = [i for i, mask in enumerate(gap_mask) if not mask]
    for i in missing_indices:
        if indices[1][i] == -1:
            trimmed_query_sequence[indices[0][i]] = "-"
    trimmed_query_sequence = "".join(trimmed_query_sequence)

    return residue_numbers, trimmed_query_sequence


def _extract_pdb_residues(
    structure: Structure, query_sequence: str) -> tuple[list[int], str, str]:
    """Extracts the sequence from the PDB file, and possibly trims the query sequence.

    Args:
        structure (PDB.PDB.Structure.Structure): The structure object from the PDB file.
        query_sequence (str): The input query sequence.

    Returns:
        tuple: The residue numbers, the trimmed query sequence, and the chain ID.
    """
    ppb = PDB.PPBuilder()
    residue_numbers = []
    trimmed_query_sequence = query_sequence

    best_match = {"score": 0, "chain": None, "chain_sequence": None}
    query_length = len(query_sequence)
    aligner = PairwiseAligner(match_score=1, mismatch_score=-2, gap_score=-1, mode="local")

    for model in structure:
        for chain in model:
            chain_sequence = "".join([str(pp.get_sequence()) for pp in ppb.build_peptides(chain)])
            chain_id = chain.id

            if query_sequence in chain_sequence:
                start_idx = chain_sequence.find(query_sequence)
                end_idx = start_idx + len(query_sequence)
                residue_numbers = [residue.id[1] for residue in list(chain)[start_idx:end_idx]]
                print(f"Found exact local correspondence between query and PDB chain {chain_id}.")
                return residue_numbers, query_sequence, chain_id

            elif chain_sequence in query_sequence:
                residue_numbers = [residue.id[1] for residue in chain]
                print(f"Found exact local correspondence between query and PDB chain {chain_id}.")
                return residue_numbers, chain_sequence, chain_id

            else:
                # compute local alignment score
                score = aligner.score(query_sequence, chain_sequence)
                score /= min(query_length, len(chain_sequence))
                if score > best_match["score"] and score > _SCORE_THRESHOLD:
                    best_match = {"score": score, "chain": chain, "chain_sequence": chain_sequence}

    if best_match["chain"] is None:
        raise ValueError(
            "Similarity between query sequence and all PDB chains below acceptance level, even"
            " with non-strict search. Please check the input files."
        )
    
    chain_id = best_match["chain"].id
    warnings.warn(
        f"Found only partial correspondence between query and PDB chain {chain_id}. "
        f"(alignment score = {best_match['score']:.2f}).\n"
        "Resulting GEMME/ESCOTT predictions may be less accurate.",
    )

    residue_numbers, trimmed_query_sequence = _process_best_match(
        best_match, query_sequence, aligner
    )

    return residue_numbers, trimmed_query_sequence, chain_id


def prepare_pdb(pdb_file: Path, alignment_file: Path) -> None:
    """Extracts sequence from PDB file and save a new PDB file with only
    the residues of the input sequence.

    Args:
        pdb_file (str): Path to the input PDB file.
        alignment_file (str): Path to the alignment file in FASTA format,
            with the query sequence as first entry.

    Returns:
        None
    """
    # Extract alignment and query sequence
    alignment = _extract_alignment(alignment_file)
    query_sequence = str(alignment[0].seq)

    # Load the PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    # Extract residues from the PDB file. If no strict match is found, try partial match
    residue_numbers, trimmed_sequence, chain_id = _extract_pdb_residues(structure, query_sequence)

    if len(trimmed_sequence.replace("-", "")) != len(residue_numbers):
        raise ValueError(
            "Trimmed query sequence length does not match PDB sequence length. Possible alignment"
            " issue."
        )

    # Save the processed PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)

    output_pdb_file = pdb_file.with_name(f"{pdb_file.stem}_processed.pdb").as_posix()
    io.save(output_pdb_file, ResidueSelect(residue_numbers, chain_id))
    print(f"Processed PDB file saved to {output_pdb_file}")

    # Trim alignment if necessary
    if trimmed_sequence != query_sequence:
        alignment = _trim_alignment(alignment, trimmed_sequence)

        # save trimmed alignment
        output_alignment_file = alignment_file.with_name(f"{alignment_file.stem}_trimmed.fasta")
        SeqIO.write(alignment, output_alignment_file, "fasta-2line")
        print(f"Trimmed alignment saved to {output_alignment_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pdb_file", type=Path, help="Path to the input PDB file.")
    parser.add_argument(
        "alignment_file",
        type=Path,
        help="""Path to the alignment file in FASTA format,
with the query sequence as first entry.
""",
    )
    args = parser.parse_args()

    prepare_pdb(args.pdb_file, args.alignment_file)
