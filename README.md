# gemme-utils

## Overview

gemme-utils is a collection of utilities for handling and processing data related to [GEMME](http://www.lcqb.upmc.fr/GEMME) software and its updated version [ESCOTT](http://gitlab.lcqb.upmc.fr/tekpinar/PRESCOTT/). These utilities include scripts for preparing input FASTA and PDB files, running Docker containers, and plotting scores on protein structures using PyMOL.


## Requirements

- Python 3.7+
- BioPython (for PDB processing)
- Docker (for running GEMME/ESCOTT)
- Pandas (for processing prediction files)
- PyMOL (for plotting scores on protein structures)


## Installation

After cloning the repository, install the required Python packages in a dedicated `conda` environment, running the following:

```bash
conda create -n gemme-env python=3.10
conda install biopython pandas
```

Alternatively, create a virtual environment:
```bash
python -m venv gemme-venv
source gemme-venv/bin/activate
pip install biopython pandas
```

To execute the scripts from any directory, source `gemme-utils.sh` file:
```bash
. /path/to/this/directory/gemme-utils.sh
```

To ensure the script loads every time you open a terminal, add the above line at the end it to the `~/.bashrc` or `~/.zshrc` file.

The command defines a new variable `GU_PATH` pointing to this directory, and a function `gemme-utils` for scripts execution. It also enables bash completion to list available scripts. Custom python scripts saved in this directory will be visible and executable as well.

Now you can call the scripts from any location with the following syntax:

```bash
gemme-utils <script.py> <arg1> ...

# Example:
cd test_data && gemme-utils prepare_alignment.py trxf2.a3m
```

Before the execution, ensure that the proper python environment is active.

### PyMOL

Official PyMOL software can be installed from [https://www.pymol.org/](https://www.pymol.org/) and it requires a valid license to be used. The provided bundles also include a version of Python.  
An open-source version of PyMOL is available at [https://github.com/schrodinger/pymol-open-source](https://github.com/schrodinger/pymol-open-source).

### Docker

Docker can be installed through the Docker Desktop application from [https://www.docker.com/](https://www.docker.com/), and it is available for Linux, Windows and MacOS.  

For some Linux distributions, it is also possible to install only Docker Engine, that allows to communicate with Docker server via command line interface (CLI). Download the appropriate version of Docker Engine from [https://docs.docker.com/engine/](https://docs.docker.com/engine/).

For Linux distributions, check out how to install Docker without root privileges at [https://docs.docker.com/engine/security/rootless/](https://docs.docker.com/engine/security/rootless/). On Windows and MacOS, instead, permissions can modified from Docker Desktop settings.


## Scripts

### prepare_alignment.py

Prepare a FASTA-like alignment for use with GEMME/ESCOTT by:
- Removing lowercase letters from sequences.
- Converting gap-like characters (`.`) into `-`.
- Removing positions containing gaps in the query sequence.
- Shortening headers to a maximum length of 50 characters.

The processed alignment is saved in the same directory as the input file.
The output file will have the same name of the input alignment file,
with `_processed.fasta` appended.

```bash
usage: prepare_alignment.py [-h] alignment_file

positional arguments:
  alignment_file  Path to the FASTA-like alignment file.

Example: python prepare_alignment.py test_data/trxf2.a3m
```


### prepare_pdb.py

Reformat PDB file such that it can be used as input for ESCOTT,
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

Requirements: BioPython

```bash
usage: prepare_pdb.py [-h] pdb_file alignment_file

positional arguments:
  pdb_file        Path to the input PDB file.
  alignment_file  Path to the alignment file in FASTA format,
                  with the query sequence as first entry.

Example: python prepare_pdb.py test_data/blat-af2 test_data/aliBLAT.fasta
```


### docker_run.py

Run GEMME or ESCOTT in a Docker container. It is possible to specify multiple input alignment files and select the number of parallel processes to use. For each alignment file, a separate Docker container will be created.  
If the Docker image is not available locally, the script will try to pull it
from Docker Hub.

Requirements: `sudo`-less Docker installed and running

```bash
usage: docker_run.py [-h] [-n NUM_PROCESSES] software alignment_files [alignment_files ...] ...

positional arguments:
  software              Software to run, either GEMME or ESCOTT (case-insensitive).
  alignment_files       Paths to the alignment files in valid FASTA format.
  additional_args       Additional arguments to pass to the Docker container 
                        in the form --key value or -k value.

options:
  -n NUM_PROCESSES, --num_processes NUM_PROCESSES
                        Number of parallel processes to use. Default: 1.
                        Note: each process will run on a separate Docker 
                        container, and will spawn additional threads during 
                        the execution. To avoid overloading the system, it 
                        is recommended to set the maximum value to the number 
                        of available CPU cores / 8.

Example: python docker_run.py escott test_data/aliBLAT.fasta -n 1 \ 
  --pdbfile test_data/blat-af2.fasta -m test_data/BLAT_ECOLX_Stiffler_2015_experimental.dat
```


### plot_on_structure.py

Plot GEMME/ESCOTT average per-position scores on the structure of the protein in a PyMOL session. The script will create a PyMOL session and it will color the structure 
based on the average per-position scores from the prediction file. 
Red color indicates positions more sensible to mutations, while blue 
indicates less sensible positions.

Requirements: PyMOL software, Pandas

```bash
usage: plot_on_structure.py [-h] [--save] prediction_file pdb_file

positional arguments:
  prediction_file  Path to the prediction file.
  pdb_file         Path to the PDB file.

options:
  -h, --help       show this help message and exit
  --save           Save the colored structure as a PSE file.

Example: python plot_on_structure.py test_data/BLAT_normPred_evolCombi_test.txt test_data/blat-af2.pdb
```


## Contact Information

For any questions, issues, or contributions, please contact gianluca.lombardi@sorbonne-universite.fr