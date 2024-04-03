import os
import time
import re
import ete3
import requests
from ete3 import faces, NodeStyle

CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "datafiles")
PATH_TO_MUSCLE_DATA = os.path.join(CWD, "datafiles", "muscle_data")
PATH_TO_UNIPROT_ENTRIES = os.path.join(CWD, "datafiles", "uniprot_entries")

FIGURE_1_RESULTS_FILENAME = "user_tree.png"
FIGURE_1_RESULTS_DIRECTORY_FILEPATH = os.path.join(os.getcwd(), "figures", "figure_1", "figure_1_results")


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

    file_path = os.path.join(FIGURE_1_RESULTS_DIRECTORY_FILEPATH, "newick.txt")
    tree = ete3.Tree(file_path)

    # change nodes to be labeled with species names instead of protein IDs
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
    ts.orientation = 0  # writes the tree left to right
    ts.show_scale = False
    ts.layout_fn = my_layout
    ts.allow_face_overlap = True
    # may or may not work

    # Render the tree to a file (e.g., tree.png)
    tree.render(os.path.join(FIGURE_1_RESULTS_DIRECTORY_FILEPATH, FIGURE_1_RESULTS_FILENAME), tree_style=ts, w=1400,
                h=1200)


def set_styles_for_each_group(tree):
    '''
    1. Using the exact names listed in included_organisms.txt, and the taxonomic groupings you created, follow Template 1 for every group of organisms you made.
    2. If some of your groups only have 1 organism in them, you need to initialize a style as detailed in Template 1. Instead of making a list and root, you will now follow Template 2 to set just the singular node to your desired style.
    3. Run the make_tree_with_species_names() function to create the tree.
    4. You can see the tree created in figure_1/figure_1_results/user_tree.png.
    5. If you make the tree and find that some colors are missing or are covering other colors, that means that the organization of the tree does not match the actual taxonomic data. At this point you can forgo the coloring, or you can follow the steps below.
    6. Identify the organisms causing the problems. These will generally be organisms not part of the same branch of the tree as the rest of the organisms in their group (the groups that you made).
    7. Color these organisms separately rather than trying to color them based on common ancestor like in Template 1. You will want to remove the problem organism from its list and make a separate code block for it based on Template 2.
    8. There is no need to make a new style for these organisms, just use the same one as the rest of their group from Template 1.
    9. If after separating all nodes that don't share close common ancestors with the rest of their groups, you still are finding color trouble, you will need to modify your code blocks.
    10. Do your best to identify which colors/code blocks are causing the issues. Sometimes this might have to be trial and error fiddling with Templates.
    11. Once you have identified the block(s) that are the problem, you need to modify the block to follow Template 3.
    12. Keep changing blocks of code to Template 3 until the problems are fixed. It will not look as nice as Template 1, so try to do this as few times as possible
    13. If you want to change the size or format of the png, the last line of code in the make_tree_with_species_names() function has arguments for w (width) and h (height) in pixels. Feel free to change them according to your preference
    '''

    '''
    # Template 1 set groups of nodes all together to one style based on common ancestor
    my_style = NodeStyle()
    my_style["bgcolor"] = "pick a color from hexadecimal system or the accepted SVG color codes at http://etetoolkit.org/docs/latest/reference/reference_treeview.html?highlight=colors#color-names"
    my_style["fgcolor"] = "black"
    Template_list = [" your", " organism", " names", " here"] # Make sure that every organism has a space before its name. Also ensure that every name is exactly the name in included_organisms.txt
    Template_root = tree.get_common_ancestor(Template_list)
    Template_root.set_style(my_style)
    '''

    '''
    # Template 2 set style for individual nodes
    my_node = tree.search_nodes(name=" organism name")[0]
    my_node.set_style(my_style)
    '''

    '''
    # Template 3 set node styles individually for a group of nodes
    my_style = NodeStyle()
    my_style["bgcolor"] = "pick a color from hexadecimal system or the accepted SVG color codes at http://etetoolkit.org/docs/latest/reference/reference_treeview.html?highlight=colors#color-names"
    my_style["fgcolor"] = "black"
    Template_list = [" some", " organism", " names", " here"] # Make sure that every organism has a space before its name. Also ensure that every name is exactly the name in included_organisms.txt
    for node_name in Template_list:
        node = tree.search_nodes(name=node_name)[0]
        node.set_style(my_style)
    '''


# Define a layout function for the tree
def my_layout(node):
    if node.is_leaf():
        # position controls the position of the text in the tree
        faces.add_face_to_node(faces.TextFace(node.name, fsize=8, fgcolor="black", fstyle="normal"), node, column=0,
                               position="branch-right")


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
    with open(os.path.join(FIGURE_1_RESULTS_DIRECTORY_FILEPATH, "newick.txt"), 'w') as outF:
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
    make_tree_with_species_names()
