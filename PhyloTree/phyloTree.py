import ete3
import os

file_path = os.path.join(os.getcwd(), "tree.txt")

tree = ete3.Tree(file_path)
ts = ete3.TreeStyle()
ts.show_leaf_name = True
tree.render("tree.png", tree_style=ts)