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
    names = open("included_organisms.txt", "r")
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
    tree.render(os.path.join(FIGURE_1_RESULTS_DIRECTORY_FILEPATH, FIGURE_1_RESULTS_FILENAME), tree_style=ts, w=1400, h=1200)

def set_styles_for_each_group(tree):
    # NOTE: This format is for K01809 run with eukaryote_shortlist.csv

    # Vertebrates
    v_style = NodeStyle()
    v_style["bgcolor"] = "lightblue"
    v_style["fgcolor"] = "black"
    Vertebrate_list = [" Canis lupus familiaris", " Rattus norvegicus", " Danio rerio", " Mus musculus", " Heterocephalus glaber", " Poecilia reticulata"]
    Vertebrate_root = tree.get_common_ancestor(Vertebrate_list)
    Vertebrate_root.set_style(v_style)


    # Fungi
    f_style = NodeStyle()
    f_style["bgcolor"] = "tan"
    f_style["fgcolor"] = "black"
    Fungi_list = [" Saccharomyces cerevisiae", " Schizosaccharomyces pombe", " Emericella nidulans"]
    # col_node = tree.search_nodes(name=colobus)[0]
    # pap_node = tree.search_nodes(name=papio)[0]
    # col_node.set_style(c_style)
    # pap_node.set_style(c_style)
    Fungi_root = tree.get_common_ancestor(Fungi_list)
    Fungi_root.set_style(f_style)

    # Protists
    p_style = NodeStyle()
    p_style["bgcolor"] = "palegoldenrod"
    p_style["fgcolor"] = "black"
    Protist_list = [" Dictyostelium discoideum"]
    ddi_node = tree.search_nodes(name=" Dictyostelium discoideum")[0]
    ddi_node.set_style(p_style)
    # Protist_root = tree.get_common_ancestor(Protist_list)
    # Protist_root.set_style(p_style)

    # Vascular Plants
    v_style = NodeStyle()
    v_style["bgcolor"] = "#a2e57b"
    v_style["fgcolor"] = "black"
    vascular_list = [" Arabidopsis thaliana", " Zea mays", " Medicago truncatula", " Oryza sativa", " Selaginella moellendorffii"]
    vascular_root = tree.get_common_ancestor(vascular_list)
    vascular_root.set_style(v_style)

    # Invertebrates
    i_style = NodeStyle()
    i_style["bgcolor"] = "#c3a1f7"
    i_style["fgcolor"] = "black"
    #F5F5DC
    Invertebrate_list = [" Drosophila melanogaster"," Caenorhabditis elegans"]
    dme_node = tree.search_nodes(name=" Drosophila melanogaster")[0]
    dme_node.set_style(i_style)
    cel_node = tree.search_nodes(name=" Caenorhabditis elegans")[0]
    cel_node.set_style(i_style)
    cin_node = tree.search_nodes(name=" Ciona intestinalis")[0]
    cin_node.set_style(i_style)
    # Invertebrate_root = tree.get_common_ancestor(Invertebrate_list)
    # Invertebrate_root.set_style(i_style)

    # Archaeplastida
    a_style = NodeStyle()
    a_style["bgcolor"] = "#ffcc80"
    a_style["fgcolor"] = "black"
    # Arch_list = ["Chlamydomonas reinhardtii"]
    # Arch_root = tree.get_common_ancestor(Arch_list)
    # Arch_root.set_style(a_style)

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
        print("Waiting for Tree...")
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
                organisms.append(organism)
            else:
                print(f"Match failed in file {fileName}")

    with open(os.path.join(os.getcwd(), "included_organisms.txt"), 'w') as outF:
                i = 0
                for organism in organisms:
                    outF.write(str(ids[i]) + " " + organism + "\n")
                    i += 1

if __name__ == "__main__":
    get_species_names_and_uniprot_ids_from_uniprot_entries()
    # make_tree_with_species_names()
