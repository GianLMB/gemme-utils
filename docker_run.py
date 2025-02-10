"""Runs GEMME or ESCOTT in a Docker container.

It is possible to specify multiple input alignment files and select
the number of parallel processes to use. For each alignment file, a
separate Docker container will be created.
If the Docker image is not available locally, the script will try to 
pull it from Docker Hub.

Requirements:
- Docker software (installed and running)
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Literal


Software = Literal["GEMME", "ESCOTT"]
_DOCKER_IMAGE = {
    "GEMME": "elodielaine/gemme:gemme", 
    "ESCOTT": "tekpinar/prescott-docker:v1.6.0"
}


class DockerRunner:
    """Handles execution of GEMME or ESCOTT within a Docker container."""

    def __init__(self, software: str, alignment_file: str, additional_args: str = ""):
        self.software = software
        self.docker_image = self._validate_software()
        alignment_file = Path(alignment_file).absolute()  # absolute path required for Docker
        self.alignment_file = alignment_file.name
        self.input_dir = alignment_file.parent
        self.args = additional_args

        self._check_docker_access()
        self._inspect_docker_image()

        print(
            f"Running {self.software} with alignment file: {self.alignment_file} in input"
            f" directory {self.input_dir}\n"
            f"Additional arguments: {self.args if self.args != '' else None}"
        )

    def _validate_software(self) -> str:
        """Validates the software selection and returns the corresponding Docker image."""
        if self.software not in _DOCKER_IMAGE:
            raise ValueError("Invalid software. Must be either 'GEMME' or 'ESCOTT'.")
        return _DOCKER_IMAGE[self.software]

    def _check_docker_access(self) -> None:
        """Ensures that Docker is accessible."""
        try:
            subprocess.run(["docker", "info"], capture_output=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                "Docker is not accessible. Ensure it is installed and running."
            ) from e

    def _inspect_docker_image(self) -> None:
        """Ensures the required Docker image is available locally, pulling it if necessary."""
        try:
            subprocess.run(
                ["docker", "inspect", self.docker_image],
                capture_output=True,
                check=True,
                text=True,
            )
            print(f"Docker image {self.docker_image} is available locally.")
        except subprocess.CalledProcessError as e:
            error_message = e.stderr.strip()
            if "No such object" in error_message:
                print(
                    f"Image {self.docker_image} not found locally. Trying pulling from Docker"
                    " Hub..."
                )
                try:
                    subprocess.run(["docker", "pull", self.docker_image], check=True)
                except subprocess.CalledProcessError as pull_error:
                    raise RuntimeError(
                        f"Error pulling Docker image: {pull_error.stderr}"
                    ) from pull_error
            else:
                raise RuntimeError(f"Error checking Docker image: {error_message}") from e

    def _build_docker_command(self) -> list[str]:
        """Constructs the Docker command to run the software."""
        # Define the base command
        mount_option = f"type=bind,source={self.input_dir},target=/project"
        base_cmd = f"docker run -ti --rm --mount {mount_option} {self.docker_image} bash -c "

        # Define the command inside the container
        if self.software == "GEMME":
            bash_cmd = (
                f"python2.7 $GEMME_PATH/gemme.py {self.alignment_file} " 
                f"-f {self.alignment_file} -r input"
            )
        else:
            bash_cmd = f"escott {self.alignment_file}"

        # Add additional arguments
        args = [arg.split("/")[-1] for arg in self.args.split()]
        bash_cmd = " ".join([bash_cmd] + args).rstrip()

        print(f"Docker command: {bash_cmd}")
        return base_cmd + f"'cd /project && {bash_cmd}'"

    def run(self) -> None:
        """Runs the software inside a Docker container."""
        cmd = self._build_docker_command()

        # Run Docker process, stream output and capture errors
        with subprocess.Popen(
            cmd, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE, text=True
        ) as process:
            try:
                # Read and print stderr output
                for line in process.stderr:
                    print(line.strip(), flush=True)

                process.wait()

                if process.returncode != 0:
                    raise RuntimeError(
                        f"Error running {self.software}. Exit code: {process.returncode}"
                    )

                print(f"{self.software} ran successfully. Output files saved in {self.input_dir}")

            except Exception as e:
                process.terminate()
                process.wait()
                raise e


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "software",
        type=str,
        help="Software to run, either GEMME or ESCOTT (case-insensitive).",
    )
    parser.add_argument(
        "alignment_files",
        type=str,
        nargs="+",
        help="Paths to the alignment files in valid FASTA format.",
    )
    parser.add_argument(
        "-n",
        "--num_processes",
        type=int,
        default=1,
        help="""Number of parallel processes to use. Default: 1.
Note: each process will run on a separate Docker 
container, and will spawn additional threads during 
the execution. To avoid overloading the system, it 
is recommended to set the maximum value to the number 
of available CPU cores / 8.
""",
    )
    parser.add_argument(
        "additional_args",
        nargs=argparse.REMAINDER,
        help="""Additional arguments to pass to the Docker container 
in the form --key value or -k value.
""",
    )
    args = parser.parse_args()
    args.software = args.software.upper()

    # get number of available CPU cores to limit the number of concurrent processes
    num_cores = os.cpu_count()
    args.num_processes = min(args.num_processes, num_cores // 8)
    print(f"Number of available CPU cores: {num_cores}. Using {args.num_processes} processes.")

    # convert additional arguments to a string
    args.additional_args = " ".join(args.additional_args)
    return args


def _run_docker(software: Software, alignment_file: str, additional_args: str) -> None:
    """Runs GEMME or ESCOTT in a Docker container"""
    runner = DockerRunner(software, alignment_file, additional_args)
    runner.run()


def run_docker(
    software: Software, alignment_files: list[str], num_processes: int, additional_args: str
) -> None:
    """Runs GEMME or ESCOTT in a Docker container, passing the alignment file
    and additional arguments as input.
    Output files will be saved in the same directory as the alignment file.
    For the list of available arguments, see the documentation of the respective
    software.
    The function does not allow to specify different arguments for each alignment file.

    Args:
        software: Software to run. Either "GEMME" or "ESCOTT".
        alignment_files: Paths to the alignment files in valid FASTA format.
        num_processes: Number of parallel processes to use.
        **kwargs: Additional arguments to pass to the Docker container. All
            additional input files should be in the same directory as the
            alignment file.

    Returns:
        None
    """
    if len(alignment_files) < num_processes:
        raise ValueError(
            "Number of processes must be less than or equal to the number of provided alignment"
            " files."
        )

    def run_single_process(alignment_file: str) -> None:
        _run_docker(software, alignment_file, additional_args)

    if num_processes > 1:
        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            executor.map(run_single_process, alignment_files)
    else:
        for alignment_file in alignment_files:
            run_single_process(alignment_file)


if __name__ == "__main__":
    args = _parse_args()
    run_docker(**vars(args))
