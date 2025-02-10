"""Plot GEMME/ESCOTT average per-position scores on the structure
of the protein in a PyMOL session.

Requirements:
- PyMOL software
- Pandas

The script will create a PyMOL session and it will color the structure 
based on the average per-position scores from the prediction file. 
Red color indicates positions more sensible to mutations, while blue 
indicates less sensible positions.
"""

try:
    import pandas as pd
except ImportError:
    raise ImportError("Pandas is required to run this script.")

import argparse
import subprocess
import time
from pathlib import Path


class PymolSession:
    """Manages a PyMOL subprocess session."""

    def __init__(self):
        cmd = ["pymol", "-p"]  # GUI mode
        self.process = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        time.sleep(2)  # Wait for PyMOL to start

    def __call__(self, command: str) -> None:
        """Sends a command to the PyMOL session."""
        self.process.stdin.write(command + "\n")
        self.process.stdin.flush()

    def close(self) -> None:
        """Closes the PyMOL session safely."""
        self.process.stdin.close()
        self.process.wait()
        self.process.terminate()


def _load_scores(prediction_file: str) -> pd.DataFrame:
    """Load the scores from the prediction file."""
    if prediction_file.endswith(".csv"): 
        scores = pd.read_csv(prediction_file, index_col=0)
    elif prediction_file.endswith(".txt"):
        scores = pd.read_csv(prediction_file, sep=r'\s+', index_col=0).transpose()
    else:
        raise ValueError("Invalid file format. Must be either CSV or TXT.")

    scores.fillna(0, inplace=True)  # set NaN values to 0
    return scores


def plot_on_structure(
    prediction_file: str, pdb_file: str, save: bool = False
) -> None:
    """Plot the average per-position scores on the structure of the protein."""
    pdb_file = Path(pdb_file)
    mol = pdb_file.stem
    
    # Load scores and compute average per-position
    scores = _load_scores(prediction_file)
    scores = scores.mean(axis=1)
    print(scores.shape)
    
    # find minimum and maximum scores
    min_score = scores.min()
    max_score = scores.max()
    if min_score >= 0 and max_score <= 1:
        print("Scores are ranksorted. Higher scores identify more important positions.")
        cmap = "blue_white_red"
    else:
        print("Scores are not ranksorted. Lower scores identify more important positions.")
        cmap = "red_white_blue"
    
    # Start PyMOL session
    pms = PymolSession()

    try:
        pms(f"load {pdb_file}")
        pms("show cartoon")

        # Set scoes as B-factors and color by spectrum
        pms(f"alter {mol}, b=10")
        pms(
            f"for res_id, score in enumerate({list(scores)}): cmd.alter(f'resi {{res_id + 1}}',"
            " f'b={score}')"
        )
        time.sleep(2)  # Wait for the scores to be set
        pms(f"spectrum b, {cmap}, minimum={min_score}, maximum={max_score}")
        
        # Add colorbar
        pms(f"ramp_new colorbar, {mol}, [{min_score}, {max_score}], {cmap.split('_')}")
        
        # Save the session if required
        if save:
            pms(f"save {mol}_colored.pse")
            pms.process.wait()
            print(f"Colored structure saved as {mol}_colored.pse")
        return pms
    
    except Exception as e:
        pms.close()
        raise e


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("prediction_file", help="Path to the prediction file.")
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument(
        "--save", action="store_true", help="Save the colored structure as a PSE file."
    )

    args = parser.parse_args()
    pms = plot_on_structure(args.prediction_file, args.pdb_file, args.save)

    try:
        print("PyMOL session is running. Press Ctrl+C to close.")
        pms.process.wait()  # Keep PyMOL open
    except KeyboardInterrupt:
        print("Closing PyMOL session.")
        pms.close()
