# bio_465_structure_function
BIO 465 - Capstone Project


## Setup
- download and install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) if you don't already have it
- add conda to PATH (if necessary)
- create a new conda environement with the target by running `conda create --name <env>` (replace `<env>` with your desired environment name)
- activate the conda environment `conda activate <env>`
- install the necessary conda packages manually if needed
   - The following script should install the appropriate packages: `conda install python requests numpy pandas biopython networkx matplotlib && conda install -c etetoolkit ete3`


<!-- [post about conditional requirements files](https://stackoverflow.com/questions/29222269/is-there-a-way-to-have-a-conditional-requirements-txt-file-for-my-python-applica) -->


## To Run
The main script for running the pipeline is found in `run_pipeline.py`. Running this file can be running using the commandline arguments of the filepath of a csv file containing a list of KEGG organism codes for target organisms and a KEGG id for a target KEGG Ortholog.
   - e.g. `python structure_function_pipeline.py "./target_organisms.csv" "K03841"`

The `launch.json` file contains a configuration that should run the pipeline for Kegg Ortholog `K00937` with the organism list found in  `diverse_target_prokaryotes.csv`.


