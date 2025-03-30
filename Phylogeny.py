from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt
from Bio import Phylo
from ete3 import Tree, TreeStyle, NodeStyle
import numpy as np
import unittest

#Inspired by https://stackoverflow.com/questions/67527658/making-a-circular-phylogenetic-tree-in-python


class PhylogeneticTreeBuilder(unittest.TestCase):    
    """
    
    Object attributes:

    Object methods:
    __init__(self, ): 

    """
    def __init__(self, alignment):
        self.calculator = DistanceCalculator("identity")  
        self.distance_matrix = self.calculator.get_distance(alignment)
        self.constructor = DistanceTreeConstructor()
        self.nj_tree = self.constructor.nj(self.distance_matrix)

        Phylo.write(self.nj_tree, "phylo_tree.nwk", "newick")
        self.builded_tree = self.build_tree_for_ete3(self.nj_tree.root)

        self.normalize_branch_lengths(self.builded_tree)
        self.ts = self.create_tree_style()
        self.apply_layout(self.builded_tree)
        
        self.builded_tree.show(tree_style=self.ts)

    def build_tree_for_ete3(self, clade):
        """ Raises an error if there are less than 1 nodes within the graph. """
        if clade.is_terminal():
            return Tree(name=clade.name, dist=clade.branch_length or 0.1)  
        ete3_tree = Tree()
        ete3_tree.dist = clade.branch_length or 0.1  
        for child in clade.clades:
            child_tree = self.build_tree_for_ete3(child)
            ete3_tree.add_child(child_tree)
        return ete3_tree

    def normalize_branch_lengths(self, tree, scale_factor=1.5):
        """ Raises an error if there are less than 1 nodes within the graph. """
        for node in tree.traverse():
            if not node.is_root():
                node.dist = max(node.dist * scale_factor, 0.1)  

    def create_tree_style(self):
        """ Raises an error if there are less than 1 nodes within the graph. """
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "c"  
        ts.arc_start = -90  
        ts.arc_span = 360  
        ts.branch_vertical_margin = 15  
        ts.scale = 200  
        return ts

    def layout_fn(self, node):
        """ Raises an error if there are less than 1 nodes within the graph. """
        nstyle = NodeStyle()
        nstyle["size"] = 10 if node.is_leaf() else 5
        nstyle["fgcolor"] = "black"
        node.set_style(nstyle)

    def apply_layout(self, tree):
        """ Raises an error if there are less than 1 nodes within the graph. """
        for n in tree.traverse():
            self.layout_fn(n)
