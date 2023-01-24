import subprocess
import sys
import json
import time
from pathlib import Path
import os

class LXPlusVanillaSubmitter:
    """
    A class for submitting vanilla job to the lxplus cluster using HTCondor.
    It takes four inputs:
    samplejson : str : The path of the input sample json
    args_string : str : The command line arguments to be passed to the analysis script
    sample_dict : dict : A dictionary of sample names and their files
    queue : str : The queue to which the job has to be submitted, default is "longlunch"
    memory : str : The memory required for the job, default is "10GB"

    The class creates a directory structure and writes job files to submit the job to the cluster. 
    The directory structure is as follows:
    1. base_dir = ".higgs_dna_vanilla_lxplus"
    2. analysis_dir = base_dir/<sample_json_name>-<current_time>
    3. input_dir = analysis_dir/inputs
    4. jobs_dir = analysis_dir/jobs

    The class has a submit method that submits all the jobs in the jobs_dir to the cluster using condor_submit.
    """
    def __init__(self, samplejson, args_string, sample_dict, queue="longlunch", memory="10GB"):
        self.original_samplejson = samplejson
        self.sample_json_name = samplejson.split("/")[-1].split(".")[0]
        self.args_string = args_string
        self.sample_dict = sample_dict
        self.queue = queue
        self.memory = memory
        current_time = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
        self.base_dir = ".higgs_dna_vanilla_lxplus"
        self.analysis_dir = os.path.join(self.base_dir, f"{self.sample_json_name}-{current_time}")
        Path(self.analysis_dir).mkdir(parents=True, exist_ok=True)
        self.input_dir = os.path.join(self.analysis_dir, "inputs")
        Path(self.input_dir).mkdir(parents=True, exist_ok=True)

        # split sample_dict in different JSON files
        self.json_files = []
        for sample in sample_dict:
            for fl in sample_dict[sample]:
                to_dump = {}
                to_dump[sample] = [fl]
                simple_file_name = fl.split("/")[-1].split(".")[0]
                json_file_name = os.path.join(self.input_dir, f"{sample}-{simple_file_name}.json")
                with open(json_file_name, "w") as jf:
                    json.dump(to_dump, jf, indent=4)
                self.json_files.append(json_file_name)

        # write job files
        self.jobs_dir = os.path.join(self.analysis_dir, "jobs")
        Path(self.jobs_dir).mkdir(parents=True, exist_ok=True)
        self.job_files = []
        current_dir = os.getcwd()
        for json_file in self.json_files:
            base_name = json_file.split("/")[-1].split(".")[0]
            job_file_name = os.path.join(self.jobs_dir, f"{base_name}.sub")
            job_file_out = os.path.join(self.jobs_dir, f"{base_name}.out")
            job_file_err = os.path.join(self.jobs_dir, f"{base_name}.err")
            with open(job_file_name, "w") as submit_file:
                arguments = self.args_string.replace(self.original_samplejson, os.path.join(current_dir, json_file)).replace(" vanilla_lxplus", " iterative")
                submit_file.write("executable = /usr/bin/env\n")
                submit_file.write(f"arguments = {sys.prefix}/bin/run_analysis.py {arguments}\n")
                submit_file.write(f"output = {job_file_out}\n")
                submit_file.write(f"error = {job_file_err}\n")
                submit_file.write(f"+RequestMemory={self.memory}\n")
                submit_file.write(f'+JobFlavour = "{self.queue}"\n')
                submit_file.write("queue 1\n")
            self.job_files.append(job_file_name)


    def submit(self):
        """
        A method to submit all the jobs in the jobs_dir to the cluster
        """
        for jf in self.job_files:
            subprocess.run(["condor_submit", jf])
        return None
