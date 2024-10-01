    
from datetime import date, datetime, timedelta
import time as time_module  # Import with a different name to avoid conflicts
import os
import sys
import glob
import multiprocessing
import subprocess

input_dir = {
    #Input and output directories 
    #Path to all scripts 
    "pathExperimentdir": "/projects/b1042/GoyalLab/aleona/20240925_BulkRNA/240925_VH01990_3_2223NG7NX/",
    #The path to the sample sheet
    "sample_sheet":"/projects/b1042/GoyalLab/aleona/20240925_BulkRNA/240925_VH01990_3_2223NG7NX/Data/Intensities/BaseCalls/SampleSheet_bcl2fastqv5.csv",
    #The path to the output directory
    "outputdir":"/projects/b1042/GoyalLab/aleona/STARAlignmentforBulk/Fastq/", #created or taken from ENSEMBLE/GENECODE 
    #Path to jobs and logs 
    "joblogpath": "/projects/b1042/GoyalLab/aleona/20240925_BulkRNA/"
}

bash_script ={
    "quest" : "True", #Default is TRUE
 ##For creating the bash scripts in Quest
    "account":"b1042",
    "partition": "genomics", 
    "nodes": "1",
    "ntasks":"8",
    "memory":"10GB",
    "time":"2:00:00", 
    "email":"AureliaLeona2028@u.northwestern.edu",
    "mail_type":"ALL" #BEGIN, END, FAIL AND ALL
}
    
    
#####################################################################################################################################
#Nothing needs to be modified in the entire pipeline beyond this.
#####################################################################################################################################

# Parse through the input files 
pathExperimentdir= input_dir["pathExperimentdir"]
sample_sheet= input_dir["sample_sheet"]
outputdir= input_dir["outputdir"] 

# Parse through the bash script inputs 
quest = bash_script['quest']
account = bash_script['account']
partition = bash_script['partition']
nodes = bash_script['nodes']
ntasks = bash_script['ntasks']
time = bash_script['time']
memory = bash_script['memory']
email = bash_script['email']
mail_type = bash_script['mail_type']

#Indicate today's date for processes
today = date.today().strftime('%Y-%m-%d')
current_time = datetime.now().strftime("%H.%M")

#Create job_directory and logs_directory 
joblogpath = input_dir["joblogpath"]
job_directory = os.path.join(joblogpath,"jobs")
os.makedirs(job_directory, exist_ok=True)
print(job_directory)
logs_directory = os.path.join(joblogpath,"logs", today)
os.makedirs(logs_directory, exist_ok=True)
print(logs_directory)

time_module.sleep(0.2)
job_file = os.path.join(job_directory, f"runbcl2fastq_{today}_{current_time}.sh")
print(job_file)
with open(job_file, 'w') as fh:
        # SLURM directives
        fh.write("#!/bin/bash\n")
        slurm_directives = [
            f"#SBATCH --account={account}",
            f"#SBATCH --nodes={nodes}",
            f"#SBATCH --partition={partition}",
            f"#SBATCH --ntasks={ntasks}",
            f"#SBATCH --time={time}",
            f"#SBATCH --mem={memory}",
            f"#SBATCH --mail-user={email}",
            f"#SBATCH --mail-type={mail_type}",
            "#SBATCH --job-name=bcl2fastq_run",
            f"#SBATCH --output={logs_directory}/%x_%j.out",
            f"#SBATCH --error={logs_directory}/%x_%j.err",
        ]
        fh.write('\n'.join(slurm_directives) + '\n\n')

        # Module loading
        fh.write("module purge\n")
        fh.write("module load bcl2fastq\n\n")

        # Path variables
        fh.write(f"SAMPLE_SHEET_PATH={sample_sheet}\n")
        fh.write(f"OUTPUT_DIR={outputdir}\n")
        fh.write(f"FILEPATH={pathExperimentdir}\n\n")

        # Change to the experiment directory
        fh.write(f"cd ${{FILEPATH}}\n\n")

        # bcl2fastq command
        bcl2fastq_cmd = [
            "bcl2fastq \\",
            "    --create-fastq-for-index-reads \\",
            "    --ignore-missing-positions \\",
            "    --ignore-missing-controls \\",
            "    --ignore-missing-filter \\",
            "    --ignore-missing-bcls \\",
            "    --output-dir=${OUTPUT_DIR} \\",
            "    --runfolder-dir=${FILEPATH} \\",
            "    --sample-sheet=${SAMPLE_SHEET_PATH}"
        ]
        fh.write('\n'.join(bcl2fastq_cmd) + '\n')

os.system("sbatch %s" %job_file)

result = subprocess.run(["sbatch", job_file], capture_output=True, text=True)
if result.returncode == 0:
    print("Job submitted successfully")
    print("Output:", result.stdout)
else:
    print("Job submission failed")
    print("Error:", result.stderr)
