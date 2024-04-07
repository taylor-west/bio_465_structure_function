import os
import time
import re
import ete3
import requests
from ete3 import faces, NodeStyle

# CWD should be figure_1 folder
CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "../../datafiles")
PATH_TO_MUSCLE_DATA = os.path.join(PATH_TO_DATAFILES, "muscle_data")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")

FIGURE_1_RESULTS_FILENAME = "eukaryote_tree.png"
FIGURE_1_RESULTS_DIRECTORY_FILEPATH = os.path.join(CWD, "figure_1_results")


def make_tree_with_species_names():
    idDict = {}
    names = open(os.path.join(PATH_TO_MUSCLE_DATA, "organisms.txt"), "r")
    lines = names.readlines()
    for line in lines:
        line = line.strip().split()
        id = line[0]
        organism = " ".join(line[1:])
        idDict[id] = str(" " + organism)
    names.close()

    get_newick_tree()

    file_path = "eukaryote_newick.txt"
    tree = ete3.Tree(file_path)

    #change nodes to be labeled with species names instead of protein IDs
    for leaf in tree.iter_leaves():
        leaf.name = idDict[leaf.name]

    node_style = NodeStyle()
    node_style["fgcolor"] = "black"
    for node in tree.traverse():
        node.set_style(node_style)

    # Traverse the tree and apply the default background color function among other style choices
    set_styles_for_each_group(tree)
    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.force_topology = False
    ts.branch_vertical_margin = 10
    ts.orientation = 0 #writes the tree left to right
    ts.show_scale = False
    ts.layout_fn = my_layout
    ts.allow_face_overlap = True
    # may or may not work

    # Render the tree to a file (e.g., tree.png)
    tree.render(os.path.join(FIGURE_1_RESULTS_DIRECTORY_FILEPATH, FIGURE_1_RESULTS_FILENAME), tree_style=ts, w=700, h=1000)

# Ustilago maydis - basidiomycota; ustilaginomycetes; ustilaginales; ustilaginaceae
# Nothobranchius furzeri - chordata; actinopteri; zyprinodontisor; noprobrachiidae;
# Trichoplax adhaerens - Placozoa; uniplacotomia; trichoplacita; trichoplacidae;
# zea mays - streptophyta; magnoliopsida; poales; poaceae;
# Branchiostoma floridae  - chordata; leptocardii; amphioxyformes; bronchiostomatidae
# Ciona intestinalis - chordata; ascidiacea; phlebobronchia; cionidae
# Macaca mulatta - chordata; mammalia; primates; cercopithecidae
# Anolis carolinensis - chordata; lepidosauria; squamata; dactyloidae
# Cryptococcus neoformans var - basidiomycota; tremellomycetes; tremellales; cryptococcaceae 2x
# Arabidopsis thaliana - streptophyta; magnoliopsida; brassicales; brassicaceae 2x
# Caenorhabditis elegans - nematoda; chromadorea; rhabditida; rhabditidae; 2x
# Dictyostelium discoideum - evosea; eumycetozoa; dictyosteliales; dictyosteliaceae
# Emericella nidulans - ascomycota; eurotiomycetes; eurotiales; aspergillaceae
# Mus musculus - chordata; mammalia; rodentia; muridae
# Neurospora crassa - ascomycota; sordariomycetes; sordariales; sordariaceae 2x
# Rattus norvegicus - chordata; mammalia; rodentia; murdiae
# Schizosaccharomyces pombe - ascomycota; schizosaccharomycetes; schizosaccharomycetales; schizosaccharomycetaceae
# Saccharomyces cerevisiae - ascomycota; saccharomycetes; saccharomycetales; saccharomycetaceae
# Danio rerio - chordata; actinopteri; cypriniformes; danionidae
# Drosophila melanogaster - arthropoda; insecta; diptera; drosophilidae

def set_styles_for_each_group(tree):
    # NOTE: This format is for K01809 run with target_eukaryotes.csv

    # Chordata
    c_style = NodeStyle()
    c_style["bgcolor"] = "#FBE665"
    c_style["fgcolor"] = "black"
    Chordata_list = [" Nothobranchius furzeri", " Branchiostoma floridae", " Macaca mulatta", " Anolis carolinensis",
                       " Mus musculus", " Rattus norvegicus", " Danio rerio"]
    Chordata_root = tree.get_common_ancestor(Chordata_list)
    Chordata_root.set_style(c_style)

    # ciona
    cio_node = tree.search_nodes(name=" Ciona intestinalis")[0]
    cio_node.set_style(c_style)


    # basidiomycota
    b_style = NodeStyle()
    b_style["bgcolor"] = "#F8D9FA"
    b_style["fgcolor"] = "black"
    b_list = [" Ustilago maydis", " Cryptococcus neoformans var", " Cryptococcus neoformans var 2"]
    Fungi_root = tree.get_common_ancestor(b_list)
    Fungi_root.set_style(b_style)

    # Placozoa
    p_style = NodeStyle()
    p_style["bgcolor"] = "#BAEDE6"
    p_style["fgcolor"] = "black"
    ddi_node = tree.search_nodes(name=" Trichoplax adhaerens")[0]
    ddi_node.set_style(p_style)

    # streptophyta
    s_style = NodeStyle()
    s_style["bgcolor"] = "lightpink" # "#FCFBB8"
    s_style["fgcolor"] = "black"
    s_list = [" Zea mays", " Arabidopsis thaliana", " Arabidopsis thaliana 2"]
    s_root = tree.get_common_ancestor(s_list)
    s_root.set_style(s_style)

    # Nematoda
    n_style = NodeStyle()
    n_style["bgcolor"] = "#E3CCA9"
    n_style["fgcolor"] = "black"
    n_list = [" Caenorhabditis elegans", " Caenorhabditis elegans 2"]
    for node_name in n_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(n_style)

    # evosea
    e_style = NodeStyle()
    e_style["bgcolor"] = "#BAF1B4"
    e_style["fgcolor"] = "black"
    e_node = tree.search_nodes(name=" Dictyostelium discoideum")[0]
    e_node.set_style(e_style)

    #ascomycota
    a_style = NodeStyle()
    a_style["bgcolor"] = "#FDCE67"
    a_style["fgcolor"] = "black"
    a_list = [" Emericella nidulans", " Neurospora crassa", " Schizosaccharomyces pombe", " Saccharomyces cerevisiae"]
    a_root = tree.get_common_ancestor(a_list)
    a_root.set_style(a_style)

    # Neurospora crassa 2
    n2_node = tree.search_nodes(name=" Neurospora crassa 2")[0]
    n2_node.set_style(a_style)

    #arthopoda
    arth_style = NodeStyle()
    arth_style["bgcolor"] = "#FCEDA8"
    arth_style["fgcolor"] = "black"
    arth_node = tree.search_nodes(name=" Drosophila melanogaster")[0]
    arth_node.set_style(arth_style)


    # # None
    # no_style = NodeStyle()
    # no_style["bgcolor"] = "white"
    # Reject_list = ["Carlito syrichta", "Microcebus murinus"]
    # for node_name in Reject_list:
    #     node = tree.search_nodes(name=node_name)[0]
    #     node.set_style(no_style)

# Define a layout function for the tree
def my_layout(node):
    if node.is_leaf():
        # position controls the position of the text in the tree
        faces.add_face_to_node(faces.TextFace(node.name, fsize=8, fgcolor="black", fstyle="normal"), node, column=0, position="branch-right")

def get_newick_tree():
    email = "alexwalbom@gmail.com"
    title = "Eukaryotes"
    with open(os.path.join(PATH_TO_MUSCLE_DATA, "alignment.aln"), 'r') as inF:
        sequence = inF.read()
    clustering = "UPGMA"
    data = {
        "email": email,
        "title": title,
        # "clustering": clustering,
        "sequence": sequence
    }
    postURL = "https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/run"
    response = requests.post(postURL, data)
    jobID = response.text
    statusURL = f"https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/status/{jobID}"
    while True:
        response = requests.get(statusURL)
        if response.text == "FINISHED":
            break
        print(f"Status: {response.text}")
        time.sleep(5)

    resultType = "tree"
    getURL = f"https://www.ebi.ac.uk/Tools/services/rest/simple_phylogeny/result/{jobID}/{resultType}"
    response = requests.get(getURL)
    newick_tree = response.text
    with open("eukaryote_newick.txt", 'w') as outF:
        outF.write(newick_tree)

def get_species_names_and_uniprot_ids_from_uniprot_entries():
    subdirectoryPath = PATH_TO_UNIPROT_ENTRIES
    fileList = os.listdir(subdirectoryPath)
    organisms = []
    ids = []
    name_dictionary = {}
    pattern = re.compile(r'OS\s+((?:(?!subsp\.)[^(.])+)?\s*')
    for fileName in fileList:
        filePath = os.path.join(subdirectoryPath, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            lines = inF.readlines()
            first_line = lines[0].split()
            ids.append(first_line[1])
            content = "".join(lines)
            match = pattern.search(content)
            if match:
                organism = match.group(1)
                organism = organism.strip()
                if organism not in name_dictionary:
                    name_dictionary[organism] = 1
                    organisms.append(organism)
                else:
                    name_dictionary[organism] += 1
                    organism += (" " + str(name_dictionary[organism]))
                    organisms.append(organism)
            else:
                print(f"Match failed in file {fileName}")

    with open(os.path.join(PATH_TO_MUSCLE_DATA, "organisms.txt"), 'w') as outF:
                i = 0
                for organism in organisms:
                    outF.write(str(ids[i]) + " " + organism + "\n")
                    i += 1

if __name__ == "__main__":
    get_species_names_and_uniprot_ids_from_uniprot_entries()
    make_tree_with_species_names()
