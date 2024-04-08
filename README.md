# BIO 465 - Capstone Project


## **Setup**



* download and install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) if you don't already have it
* add conda to PATH (if necessary)
* create a new conda environment with the target by running `conda create --name &lt;env>` (replace `&lt;env>` with your desired environment name)
* activate the conda environment `conda activate &lt;env>`
* install the necessary conda packages manually if needed
    * The following script should install the appropriate packages:  \
`conda install python requests numpy pandas biopython networkx matplotlib && conda install -c etetoolkit ete3`


## **To Run**


### **Main Pipeline**

The main script for running the pipeline is found in `run_pipeline.py`. This file can be run using the commandline arguments of the filepath of a csv file containing a list of KEGG organism codes for target organisms and a KEGG id for a target KEGG Ortholog.



* e.g. `python run_pipeline.py "K03841" "target_organisms.csv"` 

The `launch.json` file contains an IDE configuration that will run the pipeline for KEGG Ortholog `K00937` with the organism list found in `target_prokaryotes.csv`. This file can be updated to point to your desired KO ID and organism list, or you can simply run the pipeline from the command line. 


#### **Command line arguments**


##### **KEGG Ortholog ID**

The first argument that will be passed into the pipeline is the [KEGG Ortholog ID](https://www.genome.jp/kegg/ko.html) of the ortholog that you want to study. These ID's are formatted as 'K#####'. This project focused on orthologs in metabolic pathways, as these are most likely to have conserved orthologs across a broad range of organisms.

Example orthologs to use for a specific pathway can be found by:



1. Visiting [KEGG Pathways](https://www.genome.jp/kegg/pathway.html) and selecting a pathway of choice (e.g. [map01200 - 'Carbon Metabolism'](https://www.genome.jp/pathway/map01200))
2. Changing the URL to reflect a map specific to an appropriate organism. (e.g. '[https://www.genome.jp/pathway/map01200](https://www.genome.jp/pathway/map01200)' --> [https://www.genome.jp/pathway/eco01200](https://www.genome.jp/pathway/eco01200) for E. Coli)
    * (Note that the pipeline will work best when the organism that you select is found in the list of target organisms that is provided by the second command-line argument.)
3. Selecting a green highlighted arrow or node, representing a protein (e.g. [glucokinase](https://www.genome.jp/entry/eco:b2388))
4. The KO value for that protein is listed in the 'Name' box near the top of the information table (e.g. [K00845](https://www.genome.jp/entry/K00845) for glucokinase)


##### **Organism List**

The second argument passed into the pipeline is a CSV file containing information about which organisms should be used to pull proteins from when evaluating orthologs and conserved residues. The only requirement of the CSV file is a column titled `kegg_organism_code`, which contains the 3-letter code that KEGG uses to identify organisms (e.g. `eco` for E. Coli).

The example `target_prokaryotes.csv` provided with the project is a list of 50 prokaryotes that have a significant amount of phylogenetic diversity. Similarly, `target_eukaryotes.csv `is a list of 60 eukaryotes that span a significant portion of the evolutionary tree.


### **Pieces of the Pipeline**


#### **Ortholog Finder**


##### Overview

The Ortholog Finder portion of the pipeline retrieves the UniProt ID’s of the orthologs associated with any organism on the provided list. You specify which organisms you are interested in, and which ortholog you would like to target. Using this, `ortholog_identification/kegg_orthologs.py` queries the KEGG API to retrieve a list of all genes associated with the given KEGG Ortholog ID. These genes are then filtered down to only contain organisms permitted by the organism list. For each of these genes, the UniProt API is queried to convert the KEGG gene identifier to the appropriate UniProt identifier.


##### Temp Files

In the process of generating the list of UniProt ID’s, there exists the possibility that the KEGG or UniProt API could time out or experience errors. To aid in mitigating slowdowns caused by this, a reference file is created as the ID’s are processed.This file is saved to (`ortholog_identification/out/temp/temp_uniprot_ids_results.csv`). In the case that the API requests fail, a subsequent request with the same parameters will use the previous results to only make the requests that haven’t been made previously. Once all of the genes are processed correctly, the main output file (`datafiles/ortholog_uniprots/`) is created and the temporary file is destroyed.

##### Output

The main file (`ortholog_identification/kegg_orthologs.py`) returns a list of UniProt ID’s to the main pipeline (`run_pipeline.py`). It also saves a CSV file containing some key information about the generated dataframe to (`datafiles/ortholog_uniprots/`). 


##### Filtering

The `find_uniprot_gene_collections` function serves to filter out orthologs that do not have a sufficient number of highly-annotated UniProt entries. The default parameters are 10 orthologs with a minimum annotation score of 3.0.


#### **Sequence Alignment**

This portion is in the `multiple_sequence_alignment/muscle.py` file. This file takes a list of uniprot IDs generated by the ortholog finder, and retrieves each ID's entry in Uniprot as a flat-file. Those are stored in `datafiles/uniprot_entries`. It then concatenates all the flat-files together into one long string, and sends that string to the MUSCLE program via API. MUSCLE then performs a clustal alignment on the sequences of every organism we had a uniprot entry for. It returns the alignment which is then stored in `datafiles/muscle_data/alignment.aln`.

This portion of the pipeline then determines the position of every invariant residue returned in the alignment as well as the position of every invariant residue in each protein sequence. It sends that data in a dataframe to the next portion of our pipeline for further calculations and analysis.


#### **3D-Clustering**

This is in the `cluster_analysis_3d/CalculateResidueDistanceWithDataframeInput.py `file. The functions handle pulling PDB files from alphafold to match uniprot IDs, calculating average coordinates of atoms in residues, converting UniProt data to positions, extracting relevant information about active and binding sites, reading UniProt files to identify active sites and binding sites, and analyzing protein structures to find clusters of residues. The input data used was taken from the multiple sequence alignment results generated by `multiple_sequence_alignment/muscle.py` to get the expected positions of residues within the protein structure. 



* `make_expected_cluster_lists_and_find_actual_clusters() `creates a file `datafiles/cluster_data/clusters_{target_ko_id}.csv `containing a dataframe that has all the information about the clusters and associated protein sequences.

The `get_cluster_dataframe() `executes the pipeline in the script, calling other functions as needed. 


### **Figures**


#### **Figure 1**

Make sure you have the ete3 python package installed from the setup code at the top.

To replicate our Phylogenetic Tree Figures follow these steps:

Prokaryote Figure
1. Run the pipeline with `K00937 target_prokaryotes.csv` as commandline arguments
2. Go to `figures/figure_1/phylo_tree_prokaryotes.py`. 
3. Run the `main` function at the bottom of the file.
4. Go to `figures/figure_1/figure_1_results`
5. In this directory you will see the tree you have created as `prokaryote_tree.png`.

Eukaryote Figure
1. Run the pipeline with `K01809 target_eukaryotes.csv` as commandline arguments
2. Go to `figures/figure_1/phylo_tree_eukaryotes.py`. This file generates the eukaryote phylogenetic tree.
3. Run the `main` function at the bottom of the file.
4. Go to `figures/figure_1/figure_1_results`
5. In this directory you will see the tree you have created as `eukaryote_tree.png`.

**To generate your own Figure 1 Phylogenetic Tree follow these steps:**

NOTE: This will take quite a long time (at least 1.5 hours) . Please do not do this during the Reproducibility challenge. Feel free to read it though.



1. Run the pipeline with desired Kegg ID and organisms.
2. Go to `figures/figure_1/phylo_tree.py.`
3. Run only the function `get_species_names_and_uniprot_ids_from_uniprot_entries() `by commenting out `make_tree_with_species_names() `in `main() `(main is at the bottom of the file). (This takes as input the uniprot flat files in datafiles/uniprot_entries which gets overwritten each time the pipeline is run, so no arguments need to be specified to run the program. Just make sure it is run immediately after running the pipeline with the desired Kegg ID.)
4. A file should have appeared in the `figure_1` directory called `included_organisms.txt`. This file shows each organism that had enough data to be used in the pipeline.
5. The first column will be the name of the uniprot file in `datafiles/uniprot_entries` corresponding to that organism. The second column will be separated from the first by just whitespace, and it has the species names.
6. The species name for the purpose of this code normally is two words, the genus and species. However, some names in the uniprot entries do not match convention. This file uses a regex string to identify what the species name should be. To know what species name you should keep track of, we have explained the regex code below in case the need to modify it arises
    * `OS\s+((?:(?!subsp\.)[^(.])+)?\s*`
    * The function `get_species_names_and_uniprot_ids_from_uniprot_entries()` in `figures/figure_1/phylo_tree.py` includes the above regex string which becomes a pattern that matches with the species name in each uniprot entry file so the name can be extracted.
    * `OS` this part finds the line in the uniprot file that starts with OS
    * `\s+` will match any number of whitespace characters
    * `((?:(?!subsp\.)[^(.])+)?` this is the first capturing group. Everything inside will be extracted as the species name. The final `?` means it can appear zero or one time
        * `(?:(?!subsp\.)[^(.])` this is the non-capturing group. It's denoted by `?:` and is used to group together sections of text
        * `(?!subsp\.)` This is the negative lookahead denoted by `?!`. Once the search reaches `subsp.` it will terminate.
        * `[^(.]` This matches everything that is not `(` or `.`.
        * `+` means match as many characters as possible between one and infinity.
        * Altogether, this part matches as many characters as possible and will stop when it reaches `subsp.`, `(` or `.`.
    * `\s*` means match as many whitespaces as possible between zero and infinity. This is not part of the capturing group, but helps isolate it.
7. Go online to the NCBI Taxonomy Browser. [https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)
8. For every species in `included_organisms.txt` look up the taxonomy of that organism by typing the species name into the search bar
9. Keep track of the taxonomies of every organism as you go.
10. Once you have all taxonomies, decide how you want to group the organisms. Since every alignment and K0 id is different, the way you group your organisms will differ based on which had sufficient data to be included in the alignment. We grouped ours by class in prokaryote_tree.png. These groupings are for color coding purposes.
11. Once you have decided how you want to group them, open `figures/figure_1/phylo_tree.py` again.
12. Go to the `set_styles_for_each_group()` function. You will be modifying this function.
13. Follow the detailed instructions at the top of the function for how to write this code.

The Rest of Figure 1

1. The rest of our figure 1 was created with data from `datafiles/muscle_data/alignment.aln` which shows the alignment for the K0 id that was last run. 
2. The tree you generated and the data from this file were put together in Microsoft Word. The alignment was spliced and rearranged to show uniprot listed residues, match the order of the tree, and then edited to add aesthetic details. 
3. Uniprot residues used for invariant reference in alignment
    1. Eukaryotes: K01809; Mus musculus; UniProtID Q924M7; residues Q110, H112, E137, and H276
    2. Prokaryotes: K00937; Helicobacter pylori; UniProtID O25654; residues N42, Y464, R558, and H586
4. The above residues are the reference residues and all corresponding invariant residues according to the alignment are emphasized in the figure.


#### **Figure 2**

**Data used: **

_Eukaryotes_: Binding site of mannose-6-phosphate isomerase (ManA)
Homologs from: 
1. D. melanogaster (Fruit fly), UniProtID Q9VH77, shown with Q86, H88, E113, and H241 
2. S. pombe (strain 972 / ATCC 24843) (Fission yeast), UniProtID O43014, shown with Q99, H101, E126, and H265
3. D. rerio (Zebrafish) UniProtID Q3ZB95, shown with Q110, H112, E137, and H275 
4. M. musculus (Mouse), UniProtID Q924M7, shown with Q110, H112, E137, and H276

_Prokaryotes_: Binding sites of polyphosphate kinase (PPK)
Homologs from:
1. H. pylori (strain ATCC 700392 / 26695) (Campylobacter pylori), UniProtID O25654, shown with N42, Y464, R558, and H586 in the ATP binding site, as well as R372 and R401 in the Mg<sup>2+</sup> binding site
2. C. jejuni subsp. jejuni serotype O:2 (strain ATCC 700819 / NCTC 11168), UniProtID Q9PMU0, shown with N45, Y460, R553, and H580 in the ATP binding site, as well as R367 and R397 in the Mg<sup>2+</sup> binding site
3. B. anthracis, UniProtID Q81MN9, shown with N55, Y482, R578, and H606 in the ATP binding site, as well as R389 and R419 in the Mg<sup>2+</sup> binding site

**Set up Pymol:**

Install open source Pymol from GitHub: [github.com/schrodinger/pymol-open-source.git](https://github.com/schrodinger/pymol-open-source.git)

This can be done with the following conda command:

	`conda install conda-forge::pymol-open-source`

See for additional instructions if necessary:

For Windows: (https://pymolwiki.org/index.php/Windows_Install#:~:text=pymol.org/%23download-,Open%2DSource%20PyMOL,-Open%2DSource%20PyMOL)

For MacOS: (https://pymolwiki.org/index.php/MAC_Install#:~:text=PyMOL%20Users%20Archive-,Open%2DSource%20PyMOL,-Package%20managers)

Install Pymol plugin to import 3D protein structures from Alphafold2

Follow the instructions in the Readme of the following Github repository:

[github.com/APAJanssen/Alphafold2import.git](http://github.com/APAJanssen/Alphafold2import.git)

Once Pymol is installed, simply run Pymol by typing `pymol `into the command line.

Make sure to return to the Alphafold2 plugin instructions and install the plugin through the Pymol window configurations.

**Generating images with Pymol: **

(for the remainder of figure 2 creation, the commands are individually run in the Pymol window command line)

Import an individual homolog with the following command:

	`fetchAF2 [UniProt_ID]`

To get images like those included in figure 2:



1. Set the following global variables:

        ```
        set dot_density, 4
        set dot_color, gray
        ```


2. For each protein homolog, 

    	`show cartoon`, [UniProt_ID]


    	`color gray`, [UniProt_ID]


    	`select` [bind_site_name], `resi` [list of residues, concatenated by ‘+’] 


(note: if multiple proteins are loaded into the Pymol project, the above command will highlight the indicated residues in all the loaded proteins)


    	`show sticks`, [bind_site_name]


    	`show dots`, [bind_site_name]


    Example for Q9VH77:


    	`fetchAF2 Q9VH77`


            ```
            show cartoon, Q9VH77
            color gray, Q9VH77
            select Q9VH77_bind_site, resi 86+88+113+241
            show sticks, Q9VH77_bind_site
            show dots, Q9VH77_bind_site
            ```


3. On the right hand side, click the A button next to the bind_site name, then click the option labeled “center” to center the camera on the binding site.
4. Click and drag on the display window to rotate the camera until reaching a desired image.
5. Export the image by clicking File > Export image > PNG > Capture current display > Save PNG image as…


#### **Figure 3**

To generate Figure 3, run the file `figures/figure_3/make_figure_3.py` with 2 commandline arguments:



1. the target UniProt ID of the protein which you want to analyze
2. the filepath to a CSV file that contains the cluster data generated by the pipeline. This cluster data can be found in `/datafiles/cluster_data/` directory.

This script filters the cluster output data from our pipeline and generates a graphical visualization of the invariant residues and their proximity to other invariant residues in a 3D space. This graph is generated using the data for a specific protein ortholog from one organism (ie. one UniProt ID). 

It takes the cluster data produced by the 3D-Clustering section in the `cluster_analysis_3d/CalculateResidueDistanceWithDataframeInput.py` file and organizes it into communities using the NetworkX `greedy_modularity_communities` algorithm. These communities are assigned different colors in order to be more visible in the figure. These colors are assigned automatically for up to 10 communities, defined by the `COLOR_OPTIONS` constant found at the top of the `figures/figure_3/make_figure_3.py` file. \
 \
All known binding sites for the UniProt entry are then retrieved from the UniProt API and organized in a dictionary in order to compare them more easily with the community and cluster data.

Each residue/node that was predicted by the pipeline (found in the cluster data) is compared with the list of _known_ binding sites, as pulled from UniProt. If the _predicted_ (pipeline) node is found in the list of confirmed or _known_ nodes, it is assigned a label in the final graph.

A layout algorithm is chosen to decide how the nodes are visualized in the figure, and the figure is optionally annotated with additional statistics.
