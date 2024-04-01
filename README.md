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
### Main Pipeline
The main script for running the pipeline is found in `run_pipeline.py`. Running this file can be running using the commandline arguments of the filepath of a csv file containing a list of KEGG organism codes for target organisms and a KEGG id for a target KEGG Ortholog.
   - e.g. `python structure_function_pipeline.py "./target_organisms.csv" "K03841"`


#### Inputs
The `launch.json` file contains a configuration that should run the pipeline for KEGG Ortholog `K00937` with the organism list found in  `diverse_target_prokaryotes.csv`. This file can be updated to point to your desired KO ID and organism list, or you can simply run the pipeline from the command line.

   ##### KEGG Ortholog ID
   The first arguement that will be passed into the pipeline is the [KEGG Ortholog ID](https://www.genome.jp/kegg/ko.html) of the ortholog that you want to study. These ID's are formatted as 'K0#####'. This project focused on orthologs in metabolic pathways, as these are most likely to have conserved orthologs across a broad range of organisms.
   
   Example orthologs to use for a specific pathway can be found by:
   1. Visiting [KEGG Pathways](https://www.genome.jp/kegg/pathway.html) and selecting a pathway of choice (e.g. [map01200 - 'Carbon Metabolism'](https://www.genome.jp/pathway/map01200))
   2. Changing the URL to reflect a map specific to an appropriate organism. (e.g. 'https://www.genome.jp/pathway/map01200' --> https://www.genome.jp/pathway/eco01200 for E. Coli)
      - (Note that the pipeline will work best when the organism that you select is found in the list of target organisms that is provided by the sceond command-line arguement.)
   3. Selecting a green highlighted arrow or node, representing a protein (e.g. [glucokinase](https://www.genome.jp/entry/eco:b2388))
   4. The KO value for that protein is listed in the 'Name' box near the top of the information table (e.g. [K00845](https://www.genome.jp/entry/K00845) for glucokinase)

   ##### Organism List
   The second arguement passed into the pipeline is a CSV file containing information about which organisms should be considered when evaluating orthologs and conserved residues. The only requirement of the CSV file is a column titled `kegg_organism_code`, which contains the 3-letter code that KEGG uses to identify organisms (e.g. `eco` for E. Coli).

   The example `diverse_target_prokaryotes.csv` provided with the project is a list of 50 prokaryotes that have a significant amount of phylogenetic diversity.


### Figures

#### Figure 1


#### Figure 2


#### Figure 3