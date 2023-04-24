import subprocess
import sys
import json
from pathlib import Path
import os
from copy import deepcopy
import logging

logger = logging.getLogger(__name__)


class LXPlusVanillaSubmitter:
    """
    A class for submitting jobs on the CERN's LXPlus cluster using HTCondor, one job per file in a sample list of an analysis.
    The constructor creates a directory .higgs_dna_vanilla_lxplus if it does not exist and another one called .higgs_dna_vanilla_lxplus/<analysis_name>.
    If the latter exists, an exception is raised. Inside this directory two subdirectories called <inputs> and <jobs> will be created.
    In the former the split JSON files will be stored, in the latter the HTCondor related job files will be stored.

    Parameters:
        :param analysis_name: Name of the analysis.
        :type analysis_name: str
        :param analysis_dict: Dictionary containing the parameters of the analysis.
        :type analysis_dict: dict
        :param original_analysis_path: Path of the original analysis to be replaced with the new ones.
        :type original_analysis_path: str
        :param sample_dict: Dictionary containing the samples and their respective files.
        :type sample_dict: dict
        :param args_string: String containing the command line arguments.
        :type args_string: str
        :param queue: HTCondor queue to submit the job to. Defaults to "longlunch".
        :type queue: str, optional
        :param memory: Memory request for the job. Defaults to "10GB".
        :type memory: str, optional
    """

    def __init__(
        self,
        analysis_name,
        analysis_dict,
        original_analysis_path,
        sample_dict,
        args_string,
        queue="longlunch",
        memory="10GB",
    ):
        self.analysis_name = analysis_name
        self.analysis_dict = analysis_dict
        self.sample_dict = sample_dict
        self.args_string = args_string
        self.queue = queue
        self.memory = memory
        self.current_dir = os.getcwd()
        self.base_dir = os.path.join(self.current_dir, ".higgs_dna_vanilla_lxplus")
        self.analysis_dir = os.path.join(self.base_dir, self.analysis_name)
        try:
            Path(self.analysis_dir).mkdir(parents=True, exist_ok=False)
        except FileExistsError as e:
            logger.exception(
                "Directory already exists, maybe you already submitted jobs for this analysis?\n%s",
                e,
            )
            raise

        self.input_dir = os.path.join(self.analysis_dir, "inputs")
        Path(self.input_dir).mkdir(parents=True, exist_ok=True)

        # split analysis_dict and sample_dict in different JSON files
        self.json_analysis_files = []
        self.json_sample_files = []
        for sample in sample_dict:
            for fl in sample_dict[sample]:
                sample_to_dump = {}
                sample_to_dump[sample] = [fl]
                root_file_name = fl.split("/")[-1].split(".")[0]
                sample_file_name = os.path.join(
                    self.input_dir, f"{sample}-{root_file_name}.json"
                )
                with open(sample_file_name, "w") as jf:
                    json.dump(sample_to_dump, jf, indent=4)
                self.json_sample_files.append(sample_file_name)
                an_file_name = os.path.join(
                    self.input_dir, f"AN-{sample}-{root_file_name}.json"
                )
                an_to_dump = deepcopy(self.analysis_dict)
                an_to_dump["samplejson"] = sample_file_name
                with open(an_file_name, "w") as jf:
                    json.dump(an_to_dump, jf, indent=4)
                self.json_analysis_files.append(an_file_name)

        # write job files
        self.jobs_dir = os.path.join(self.analysis_dir, "jobs")
        Path(self.jobs_dir).mkdir(parents=True, exist_ok=True)
        self.job_files = []
        for json_file in self.json_analysis_files:
            base_name = json_file.split("/")[-1].split(".")[0]
            job_file_name = os.path.join(self.jobs_dir, f"{base_name}.sub")
            job_file_out = os.path.join(self.jobs_dir, f"{base_name}.out")
            job_file_err = os.path.join(self.jobs_dir, f"{base_name}.err")
            with open(job_file_name, "w") as submit_file:
                arguments = self.args_string.replace(
                    original_analysis_path, json_file
                ).replace(" vanilla_lxplus", " iterative")
                submit_file.write("executable = /usr/bin/env\n")
                submit_file.write(
                    f"arguments = {sys.prefix}/bin/run_analysis.py {arguments}\n"
                )
                submit_file.write(f"output = {job_file_out}\n")
                submit_file.write(f"error = {job_file_err}\n")
                submit_file.write(f"request_memory = {self.memory}\n")
                submit_file.write("getenv = True\n")
                submit_file.write(f'+JobFlavour = "{self.queue}"\n')
                submit_file.write("queue 1\n")
            self.job_files.append(job_file_name)

    def submit(self):
        """
        A method to submit all the jobs in the jobs_dir to the cluster
        """
        for jf in self.job_files:
            if self.current_dir.startswith("/eos"):
                # see https://batchdocs.web.cern.ch/troubleshooting/eos.html#no-eos-submission-allowed
                subprocess.run(["condor_submit", "-spool", jf])
            else:
                subprocess.run(["condor_submit", jf])
        return None
