import os
import time
import re
import ete3
import requests
from ete3 import faces, NodeStyle


CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "../../datafiles")
PATH_TO_MUSCLE_DATA = os.path.join(CWD, PATH_TO_DATAFILES, "muscle_data")
PATH_TO_UNIPROT_ENTRIES = os.path.join(CWD, PATH_TO_DATAFILES, "uniprot_entries")

FIGURE_1_RESULTS_FILENAME = "figure1.png"
FIGURE_1_RESULTS_FILE_FILEPATH = os.path.join(os.getcwd(), FIGURE_1_RESULTS_FILENAME)


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

    #get_newick_tree()

    file_path = os.path.join(os.getcwd(), "newick.txt")
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
    tree.render(os.path.join(os.getcwd(), "figure_1_results/prokaryote_tree.png"), tree_style=ts, w=700, h=1000)

# psuedomoadota / gammaproteobacteria / enterobacterales (eco) Escherichia coli
# actinomycetota / actinomycetes / mycobacteriales (mtu) Mycobacterium tuberculosis
# pseudomonadota / gammaproteobacteria / vibrionales (vch) Vibrio cholerae serotype 01
# campylobacterota / epsilonproteobacteria / campylobacterales (hpy) Helicobacter pylori
# pseudomonadota / gammaproteobacteria / pseudomonadales / (pae) Pseudomonas aeruginosa
# bacillota / bacilli / bacillales (ban) Bacillus anthracis
# pseudomonadota / gammaproteobacteria / enterobacterales (ype) Yersinia pestis
# pseudomonadota / proteobacteria / neisseriales (ngo) Neisseria gonorrhoeae
# basillota / Baccilli / bacillales (sep) Staphylococcus epidermidis
# campylobacterota / epsilonprteobacteria / campylobacterales (cje) Campylobacter jejuni
# pseudomonadota / gammaproteobacteria / enterobacterales (sfl) Shigella flexneri
# Euryarchaeota /  methanobacteria / methanobacteriales (msm) Mycolicibacterium smegmatis
# cyanobacetiota / cyanophyceae / nostocales  (ana) Nostoc sp
# aquificota / aquificae / aquificales (aeo) Aeromonas salmonicida
# pseudommonadota / gammaprotetobacteria / enterobacterales (sen) Saccharopolyspora erythraea
# peudomonadota / gammaproteobacteria / moraxellales (aba) Koribacter versatilis
# spirochaetota / spirochatia / spirochatales (tdn) Sulfurimonas denitrificans
# bacteroidota / bacteroidia / bacteroidales (bfr) Bacteroides fragilis
# pseudomonadota / betaproteobacteria / burkaholderiales (bpe) Bordetella pertussis
# pseudomondadota / gammaproteobacteria / enterobacterales (kpn) Klebsiella pneumoniae
# pseudomonadota / alphaproteobacteria / hyphomicrobiales (rle) Rhizobium johnstonii
# actinomycetota / actinomycetes / propionibacteriales (pac) Cutibacterium acnes

def set_styles_for_each_group(tree):
    # NOTE: This format is for K00937 run with diverse_target_prokaryotes.csv

    # White(  # FFFFFF)
    #     Pale
    # Blue(  # B0E0E6)
    #     Light
    # Green(  # 90EE90)
    #     Lavender(  # E6E6FA)
    #         Peach(  # FFDAB9)
    #             Pale
    # Yellow(  # FFFFE0)
    #     Misty
    # Rose(  # FFE4E1)
    #     Light
    # Gray(  # D3D3D3)

    # Bacilli
    bci_style = NodeStyle()
    bci_style["bgcolor"] = "#B0E0E6" #"lightblue"
    bci_style["fgcolor"] = "black"
    Bacilli_list = [" Bacillus anthracis", " Staphylococcus epidermidis"]
    Bacilli_root = tree.get_common_ancestor(Bacilli_list)
    Bacilli_root.set_style(bci_style)


    # gammaproteobacteria (enterobacterales)
    g_style = NodeStyle()
    g_style["bgcolor"] = "lightsteelblue" #"tan"
    g_style["fgcolor"] = "black"
    Gamma_list = [" Escherichia coli", " Yersinia pestis", " Shigella flexneri", " Klebsiella pneumoniae"]
    Gamma_root = tree.get_common_ancestor(Gamma_list)
    Gamma_root.set_style(g_style)

    # gammaproteobacteria (Other)
    # g_style = NodeStyle()
    # g_style["bgcolor"] = "tan"
    # g_style["fgcolor"] = "black"
    GammaO_list = [" Pseudomonas aeruginosa", " Koribacter versatilis"]
    for node_name in GammaO_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(g_style)

    # actinomycetes
    act_style = NodeStyle()
    act_style["bgcolor"] = "#E6E6FA" #"palegoldenrod"
    act_style["fgcolor"] = "black"
    Actinomycetes_list = [" Mycobacterium tuberculosis", " Cutibacterium acnes"]
    for node_name in Actinomycetes_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(act_style)

    # epsilonproteobacteria
    e_style = NodeStyle()
    e_style["bgcolor"] = "mediumslateblue" #"#a2e57b"
    e_style["fgcolor"] = "black"
    Epsilon_list = [" Helicobacter pylori", " Campylobacter jejuni"]
    for node_name in Epsilon_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(e_style)

    # betaproteobacteria
    b_style = NodeStyle()
    b_style["bgcolor"] = "lightcyan" #"#ffcc80"
    b_style["fgcolor"] = "black"
    Beta_list = [" Bordetella pertussis", " Neisseria gonorrhoeae"]
    Beta_root = tree.get_common_ancestor(Beta_list)
    Beta_root.set_style(b_style)

    # bacteroides
    bac_style = NodeStyle()
    bac_style["bgcolor"] = "#FFE4E1" #"#c3a1f7"
    bac_style["fgcolor"] = "black"
    Bacteroides_list = [" Bacteroides fragilis", " Bacteroides fragilis 2"]  # there are two of these so this might cause issues...
    for node_name in Bacteroides_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(bac_style)

    # None
    no_style = NodeStyle()
    no_style["bgcolor"] = "white"
    no_style["fgcolor"] = "black"
    Reject_list = [" Mycolicibacterium smegmatis", " Nostoc sp", " Aeromonas salmonicida", " Sulfurimonas denitrificans"]
    for node_name in Reject_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(no_style)

    # loner
    node_name = " Saccharopolyspora erythraea"
    node = tree.search_nodes(name=node_name)[0]
    node.set_style(g_style)

    #loner
    node_name = " Vibrio cholerae serotype O1"
    node = tree.search_nodes(name=node_name)[0]
    node.set_style(g_style)

    # alphaproteobacteria
    a_style = NodeStyle()
    a_style["bgcolor"] = "#D3D3D3" #"gold"
    a_style["fgcolor"] = "black"
    Alpha_node = " Rhizobium johnstonii"
    node = tree.search_nodes(name=Alpha_node)[0]
    node.set_style(a_style)





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
        #"clustering": clustering,
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
    with open(os.path.join(os.getcwd(), "newick.txt"), 'w') as outF:
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
                if organism not in name_dictionary:
                    name_dictionary[organism] = 1
                    organisms.append(organism)
                else:
                    name_dictionary[organism] += 1
                    organism += str(name_dictionary[organism])
                    organisms.append(organism)
            else:
                print(f"Match failed in file {fileName}")

    with open(os.path.join(PATH_TO_MUSCLE_DATA, "organisms.txt"), 'w') as outF:
                i = 0
                for organism in organisms:
                    outF.write(str(ids[i]) + " " + organism + "\n")
                    i += 1

if __name__ == "__main__":
    #get_species_names_and_uniprot_ids_from_uniprot_entries()
    make_tree_with_species_names()