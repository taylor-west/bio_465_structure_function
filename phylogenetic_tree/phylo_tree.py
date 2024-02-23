import os

import ete3
from ete3 import faces, NodeStyle


def make_tree_with_species_names():
    idDict = {}
    names = open("organisms.txt", "r")
    lines = names.readlines()
    for line in lines:
        line = line.strip().split()
        id = line[0]
        organism = " ".join(line[1:])
        idDict[id] = organism
    names.close()

    file_path = os.path.join(os.getcwd(), "tree.txt")
    tree = ete3.Tree(file_path)

    #change nodes to be labeled with species names instead of protein IDs
    for leaf in tree.iter_leaves():
        leaf.name = idDict[leaf.name]

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

    # Render the tree to a file (e.g., tree.png)
    tree.render("tree.png", tree_style=ts)

def set_styles_for_each_group(tree):
    # Hominidae
    h_style = NodeStyle()
    h_style["bgcolor"] = "lightblue"
    Hominidae_list = ["Homo sapiens", "Pan troglodytes", "Gorilla gorilla"]
    Hominidae_root = tree.get_common_ancestor(Hominidae_list)
    Hominidae_root.set_style(h_style)


    # Cercopithecidae
    c_style = NodeStyle()
    c_style["bgcolor"] = "#c3a1f7"
    Cercopithecidae_list = ["Macaca nemestrina", "Macaca mulatta"]
    colobus = "Colobus angolensis"
    papio = "Papio anubis"
    col_node = tree.search_nodes(name=colobus)[0]
    pap_node = tree.search_nodes(name=papio)[0]
    col_node.set_style(c_style)
    pap_node.set_style(c_style)
    Cercopithecidae_root = tree.get_common_ancestor(Cercopithecidae_list)
    Cercopithecidae_root.set_style(c_style)

    # Tarsiidae
    t_style = NodeStyle()
    t_style["bgcolor"] = "#FFFFFF"
    Tarsiidae_list = ["Carlito syrichta", "Microcebus murinus"]
    Tarsiidae_root = tree.get_common_ancestor(Tarsiidae_list)
    Tarsiidae_root.set_style(t_style)

    # Cebidae
    ceb_style = NodeStyle()
    ceb_style["bgcolor"] = "#a2e57b"
    Cebidae_list = ["Callithrix jacchus", "Saimiri boliviensis"]
    Cebidae_root = tree.get_common_ancestor(Cebidae_list)
    Cebidae_root.set_style(ceb_style)

    # Sciuridae
    s_style = NodeStyle()
    s_style["bgcolor"] = "#fff59d"
    Sciuridae_list = ["Heterocephalus glaber", "Castor canadensis","Dipodomys ordii"]
    Sciuridae_root = tree.get_common_ancestor(Sciuridae_list)
    Sciuridae_root.set_style(s_style)

    # Muridae
    m_style = NodeStyle()
    m_style["bgcolor"] = "#ffcc80"
    Muridae_list = ["Rattus norvegicus", "Mus musculus", "Mus caroli", "Cricetulus griseus", "Nannospalax galili", "Mesocricetus auratus"]
    Muridae_root = tree.get_common_ancestor(Muridae_list)
    Muridae_root.set_style(m_style)

    # None
    no_style = NodeStyle()
    no_style["bgcolor"] = "white"
    Reject_list = ["Carlito syrichta", "Microcebus murinus"]
    for node_name in Reject_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(no_style)

# Define a layout function for the tree
def my_layout(node):
    if node.is_leaf():
        # position controls the position of the text in the tree
        faces.add_face_to_node(faces.TextFace(node.name, fsize=8, fgcolor="black", fstyle="italic"), node, column=0, position="branch-right")

def get_species_names_and_uniprot_ids_from_uniprot_entries():
    currentDirectory = os.getcwd()
    folderName = "uniprot_entries"
    subdirectoryPath = os.path.join(currentDirectory, folderName)
    fileList = os.listdir(subdirectoryPath)
    organisms = []
    ids = []
    for fileName in fileList:
        filePath = os.path.join(subdirectoryPath, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            lines = inF.readlines()
            first_line = lines[0].split()
            ids.append(first_line[1])
            for line in lines:
                if line.startswith("OS"):
                    words = line.split()
                    organism = words[1] + " " + words[2]
                    organisms.append(organism)
                    break

    with open(os.path.join(currentDirectory, "organisms.txt"), 'w') as outF:
                i = 0
                for organism in organisms:
                    outF.write(str(ids[i]) + " " + organism + "\n")
                    i += 1

if __name__ == "__main__":
    get_species_names_and_uniprot_ids_from_uniprot_entries()
    make_tree_with_species_names()
