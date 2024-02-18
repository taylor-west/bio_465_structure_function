import os

import ete3
from ete3 import faces, NodeStyle


def make_tree_with_species_names():
    idDict = {}
    names = open("organisms.txt", "r")
    lines = names.readlines()
    for line in lines:
        line = line.strip().split()
        idDict[line[0]] = line[1]
    names.close()

    file_path = os.path.join(os.getcwd(), "tree.txt")
    tree = ete3.Tree(file_path)

    #change nodes to be labeled with species names instead of protein IDs
    for leaf in tree.iter_leaves():
        leaf.name = idDict[leaf.name]

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.force_topology = False
    # Hominidae
    h_style = NodeStyle()
    h_style["bgcolor"] = "#9ac4f8"
    Hominidae_list = ["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla"]
    Hominidae_root = tree.get_common_ancestor(Hominidae_list)
    Hominidae_root.set_style(h_style)

    # Cercopithecidae
    c_style = NodeStyle()
    c_style["bgcolor"] = "#c3a1f7"
    Cercopithecidae_list = ["Papio_anubis", "Macaca_nemestrina",
                            "Colobus_angolensis", "Macaca_mulatta"]
    Cercopithecidae_root = tree.get_common_ancestor(Cercopithecidae_list)
    Cercopithecidae_root.set_style(c_style)

    # Tarsiidae
    t_style = NodeStyle()
    t_style["bgcolor"] = "#FFFFFF"
    Tarsiidae_list = ["Carlito_syrichta", "Microcebus_murinus"]
    Tarsiidae_root = tree.get_common_ancestor(Tarsiidae_list)
    Tarsiidae_root.set_style(t_style)

    # Cebidae
    ceb_style = NodeStyle()
    ceb_style["bgcolor"] = "#a2e57b"
    Cebidae_list = ["Callithrix_jacchus", "Saimiri_boliviensis"]
    Cebidae_root = tree.get_common_ancestor(Cebidae_list)
    Cebidae_root.set_style(ceb_style)

    # Heteromyidae
    #het_style = NodeStyle()
    #het_style["bgcolor"] = "#ffb3b3"
    #Heteromyidae_list = ["Dipodomys_ordii"]
    #Heteromyidae_root = tree.get_common_ancestor(Heteromyidae_list)
    #Heteromyidae_root.set_style(het_style)

    # Sciuridae
    s_style = NodeStyle()
    s_style["bgcolor"] = "#fff59d"
    Sciuridae_list = ["Heterocephalus_glaber", "Castor_canadensis","Dipodomys_ordii"]
    Sciuridae_root = tree.get_common_ancestor(Sciuridae_list)
    Sciuridae_root.set_style(s_style)

    # Cheirogaleidae
    #Che_style = NodeStyle()
    #Che_style["bgcolor"] = "#ffcc80"
    #Cheirogaleidae_list = ["Microcebus_murinus"]
    #Cheirogaleidae_root = tree.get_common_ancestor(Cheirogaleidae_list)
    #Cheirogaleidae_root.set_style(Che_style)

    # Muridae
    m_style = NodeStyle()
    m_style["bgcolor"] = "#ffcc80"
    Muridae_list = ["Rattus_norvegicus", "Mus_musculus", "Mus_caroli", "Cricetulus_griseus", "Nannospalax_galili", "Mesocricetus_auratus"]
    Muridae_root = tree.get_common_ancestor(Muridae_list)
    Muridae_root.set_style(m_style)

    ts.layout_fn = my_layout
    # Render the tree to a file (e.g., tree.png)
    tree.show()
    tree.render("tree.png", tree_style=ts)


# Define a layout function to align leaf names to the right
def my_layout(node):
    if node.is_leaf():
        faces.add_face_to_node(faces.TextFace(node.name, fsize=10, fgcolor="black", fstyle="italic"), node, column=0, position="branch-right")

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
                    organism = words[1] + "_" + words[2]
                    organisms.append(organism)
                    break

    with open(os.path.join(currentDirectory, "organisms.txt"), 'w') as outF:
                i = 0
                for organism in organisms:
                    outF.write(str(ids[i]) +  " " + organism + "\n")
                    i += 1

if __name__ == "__main__":
    get_species_names_and_uniprot_ids_from_uniprot_entries()
    make_tree_with_species_names()
