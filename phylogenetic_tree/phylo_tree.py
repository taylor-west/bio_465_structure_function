import os
from ete3 import NodeStyle

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

def get_species_names_and_uniprot_ids_from_uniprot_entries():
    currentDirectory = os.getcwd()
    folderName = "uniprotEntries"
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
