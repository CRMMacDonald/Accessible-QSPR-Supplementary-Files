# Accessible-QSPR-Supplementary-Files
Data to accompany the paper "An accessible property classification framework to predict the solubility of functionalised naphthalenes and rylenes in organic solvents".

# Enviromment Installation
## Installation of virtual environment

Create an isolated environment to avoid conflict with system packages:

Create a working folder:
$ mkdir program_name && cd program_name

Clone the repository:
$ git clone https://github.com/CRMMacDonald/Accessible-QSPR-Supplementary-Files.git

Create a virtual environment, using either conda, python venv, or ideally providing both options, for example:
$ python -m venv prog_name_venv

Activate it:
$ source ./prog_name_venv/bin/activate

## Installation of dependencies and packages:

(prog_name_venv) $ python -m pip install -r requirements.txt

Coda option (via the environment.yml file):

(prog_name_venv) $ conda env create -f environment.yml
(prog_name_venv) $conda activate prog_name_venv
