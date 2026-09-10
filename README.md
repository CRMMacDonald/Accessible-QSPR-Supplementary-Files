# Accessible QSPR Supplementary Files
Supplementary data and notebooks accompanying the paper:

*"An accessible property classification framework to predict the solubility of functionalised naphthalenes and rylenes in organic solvents".*

## Environment Installation
### Creation of Virtual Environment

Create an isolated environment to avoid conflicts with system packages.

Create a working folder:
```bash
mkdir program_name && cd program_name
```

Clone the repository:
```bash
git clone https://github.com/CRMMacDonald/Accessible-QSPR-Supplementary-Files.git
```

Create a virtual environment using Python venv:
```bash
python -m venv venv_name
```

Activate it:
```bash
source ./venv_name/bin/activate
```

### Installation of Dependencies and Packages
With the environment activated, the required dependencies and packages can be installed using pip or Conda.

Pip option (via the requirements.txt file):
```bash
python -m pip install -r requirements.txt
```

Conda option (via the environment.yml file):

```bash
conda env create -f environment.yml
conda activate venv_name
```

## Jupyter Notebooks

The repository contains the following notebooks:

### [Smiles Writer Explanation Notebook.ipynb](https://github.com/CRMMacDonald/Accessible-QSPR-Supplementary-Files/blob/main/Smiles%20Writer%20Explanation%20Notebook.ipynb)

Using a CSV file containing SMILES representations of amino acid side chains (such as those provided in AASC_SMILES_Curated.csv) and a core SMILES structure.

The notebook outputs a CSV with each possible output SMILES structure where a single wildcard is substituted.

### [Dipeptide Smiles Writer Explanation Notebook.ipynb](https://github.com/CRMMacDonald/Accessible-QSPR-Supplementary-Files/blob/main/Dipeptide%20Smiles%20Writer%20Explanation%20Notebook.ipynb)

Using a CSV file containing SMILES representations of amino acid side chains (such as those provided in AASC_SMILES_Curated.csv) and a core SMILES structure.

The notebook outputs a CSV with each possible output SMILES structure where every possible pairwise permutation of wildcards is substituted.
