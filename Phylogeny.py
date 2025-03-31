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
    A class with the purpose of creating a phylogenetic tree of a Multiple Sequence Alignment.
    
    Object attributes:
        - alignment: A Multiple Sequence Alignment. 


    Object methods:
        - __init__(self, alignment): Creates an ete3 Tree and plots it.
        - build_tree_for_ete3(self, clade): Convert a biopython tree to an ete3 tree.
        - normalize_branch_lengths(self, tree, scale_factor=1.5):
        - create_tree_style(self):
        - layout_fn(self, node):
        - apply_layout(self, tree):
    """

    def __init__(self, alignment):

        #From the alignment, create a distance matrix
        self.calculator = DistanceCalculator("identity")  
        self.distance_matrix = self.calculator.get_distance(alignment)
        #From the distance matric, create a neighbour joining tree
        self.constructor = DistanceTreeConstructor()
        self.nj_tree = self.constructor.nj(self.distance_matrix)

        #Build and customise the ete3 Tree
        self.builded_tree = self.build_tree_for_ete3(self.nj_tree.root)
        self.normalize_branch_lengths(self.builded_tree)
        self.ts = self.create_tree_style()
        self.apply_layout(self.builded_tree)

        #Show the ete3 tree
        self.builded_tree.show(tree_style=self.ts)


    def build_tree_for_ete3(self, clade) -> Tree:
        """
        Convert a biopython tree to an ete3 tree.

        Args:
            - clade (Bio.Phylo.BaseTree.Clade): a clade from a Bio.Phylo tree structure to be converted into an ete3 tree. 

        Return:
            - ete3.Tree: an ete3 Tree translated from Bio.Phylo tree structure.
        """
        #Check if the clade is a leaf and if so, return a tree with that clade's name and defined branch length.
        #This part of the code is essential for the recursive looping of that function.
        if clade.is_terminal():
            return Tree(name=clade.name, dist=clade.branch_length or 0.1)  
        
        #Initialise an ete3 Tree and define a branch length
        ete3_tree = Tree()
        ete3_tree.dist = clade.branch_length or 0.1  

        #For the current leaf node, check for all of its children node, make a tree of them recursively, and add it as a subtree
        for child in clade.clades:
            child_tree = self.build_tree_for_ete3(child)                        #Build a tree with the children nodes recursively
            ete3_tree.add_child(child_tree)                                     #Add that tree as a subtree

        return ete3_tree


    def normalize_branch_lengths(self, tree, scale_factor=1.5) -> None:
        """
        Traverses through the tree's nodes and normalises their branch lengths

        Args:
            - tree: an ete3 tree with nodes and branches
            - scale_factor (int): a constant scale to apply to all branch lengths.
        """
        #For each node, if it is not root, normalise its branch length, but with a min of 0.1.
        for node in tree.traverse():
            if not node.is_root():
                node.dist = max(node.dist * scale_factor, 0.1)  


    def create_tree_style(self) -> TreeStyle:
        """
        Customises the ete3 tree.
        """
        ts = TreeStyle()                            #Create a TreeStyle object to customise the ete3 tree.
        ts.show_leaf_name = True                    #Show the leaf names.
        ts.mode = "c"                               #Set the tree layout to circular
        ts.arc_start = -90                          #Start the layout from -90 degrees. Which is looking up.
        ts.arc_span = 360                           #Span the tree 360 degrees
        ts.branch_vertical_margin = 15              #Set spacing between branches
        ts.scale = 200                              #Scale the branch lengths.
        return ts


    def layout_fn(self, node) -> None:
        """
        Creates a custom layout to be applied to all nodes.

        Args:
            - node: each node withing the tree.
        """
        nstyle = NodeStyle()                                    #Create a NodeStyle object to customise the ete3 tree.
        nstyle["size"] = 10 if node.is_leaf() else 5            #Set the size of leafes to be 10 and for internal nodes to be 5.
        nstyle["fgcolor"] = "black"                             #Set the color to black
        node.set_style(nstyle)                                  #Set the style to the custom parameters we have chosen


    def apply_layout(self, tree) -> None:
        """
        Aplies the layout_fn() style to the tree.

        Args:
            - tree: an ete3 tree with nodes and branches
        """
        #Traverse through each node and apply the node style.
        for n in tree.traverse():
            self.layout_fn(n)

