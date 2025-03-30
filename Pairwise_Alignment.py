from Bio import pairwise2
import numpy as np
import unittest


class PairwiseAlignment(unittest.TestCase):
    """
    A class with the purpose of making pairwise alignments between a single sequence and each sequence withing a list of sequences
    
    Object attributes:
        - target_record (str): the targeted sequence, each alignment would be done with respect to it
        - breeds_records (list[str, ...]): a list of sequences, each sequence would be aligned with the targeted sequence seperetaly. 


    Object methods:
        - __init__(self, graph): Create a PairwiseAlignment object
        - difference_percent(self, alignments): Calculates the difference in percents between the two records from a pairwise alignment
        - single_pairwise_alignment(self, method = "global", index = 0): Does a single pairwise alignment between the targeted record 
            and the best aligned sequence from breeds records.
        - all_pairwise_alignments(self): Aligns each record withing breeds_records with the target_record in order to find the best alignment 
            between them, using their scores as a measurement
        - best_index_get(self): Returns the index of the best aligned sequence within the self.breeds_records.
        - best_score_get(self): Returns the best score of the best aligned sequence between the self.target_record and self.breeds_records.
        - scores_get(self): Returns a list of all scores of alignment between the self.target_record and self.breeds_records.
    """
    def __init__(self, target_record, breeds_record):

        self.target_record = target_record                                  #The targeted record
        self.breeds_records = breeds_record                                 #The list of records to align the targeted sequence with
        self.scores = np.zeros(len(self.breeds_records))                    #The scores of alignments, assigned with 0 as placeholders.
        self.best_score = 0                                                 #The best score from all alignments
        self.best_index = 0                                                 #The location of the sequence best aligned with the targeted sequence


    def difference_percent(self, alignments) -> float:
        """ 
        Calculates the difference in percents between the two records from a pairwise alignment 
        
        Args:
            - alignments (list[tuple[str, str,float,float,float]]): the possibe pairwise alignments between two sequences
        Return:
            - float: the percentage difference between the two sequences that were aligned
        """
        differences = 0
        
        #Choose the first possible alignment from the list of all alignments
        seq1 = alignments[0][0]                             #The first sequence of the first possible alignment
        seq2 = alignments[0][1]                             #The second sequence of the first possible alignment

        #Loop through each letter for both sequences through their alignment
        for i in range(len(seq1)):
            #For any difference, add a counter
            if seq1[i] != seq2[i]:
                differences += 1

        #Calculate, and return, the difference as a percentage
        return (differences / len(alignments[0][0])) * 100
    
    
    def single_pairwise_alignment(self, method = "global", index = 0) -> list[tuple[str, str,float,float,float]]:
        """ 
        Does a single pairwise alignment between the targeted record and the best aligned sequence from breeds records.

        Args:
            - method (str): the method of alignment. It can be "local" or the default which is global alignment.
            - index (int): the location of the sequence to be aligned with the targeted sequence
        Return:
            - list[tuple[str, str,float,float,float]]: a list of possible alignments. together with their aligned sequences and their scores.
        """
        #If the chosen method is "local", use local pairwise alignment
        if method == "local":
            return pairwise2.align.localxx(self.target_record.seq, self.breeds_records[index].seq)
        #The default, use global pairwise alignment between the two chosen sequences
        return pairwise2.align.globalxx(self.target_record.seq, self.breeds_records[index].seq)


    def all_pairwise_alignments(self) -> None:
        """ 
        Aligns each record withing breeds_records with the target_record in order to find the best alignment between them, 
        using their scores as a measurement
        """
        #Loop through all sequences within the list of breeds records.
        for i, breed in enumerate(self.breeds_records):
            #Align the targeted sequence with the i-th sequence withing the breeds records
            self.scores[i] = pairwise2.align.globalxx(self.target_record.seq, breed.seq, score_only=True)

            #Either there were no other recorded scores, or a better one is found. When that happens, update the best saved aligned sequence.
            if self.best_score is None or self.best_score <= self.scores[i]:
                self.best_score = self.scores[i]
                self.best_index = i  

            #Print a process bar, for convenience.
            print(f"Processed {(i + 1)}/{len(self.breeds_records)} breeds ({((i + 1)/len(self.breeds_records))*100:.2f}%)", end="\r")


    def best_index_get(self) -> int:
        """ 
        Returns the index of the best aligned sequence within the self.breeds_records. 
        """
        return self.best_index
    

    def best_score_get(self) -> int:
        """ 
        Returns the best score of the best aligned sequence between the self.target_record and self.breeds_records. 
        """
        return self.best_score
    

    def scores_get(self) -> list:
        """ 
        Returns a list of all scores of alignment between the self.target_record and self.breeds_records. 
        """
        return self.scores