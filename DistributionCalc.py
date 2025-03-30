import numpy as np
import unittest
from scipy.stats import genextreme 


##The code bellow is inspired by the statistics of sequence similarity scores done by Blast: https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html


class DistributionCalc(unittest.TestCase):
    """
    A class with the purpose of making pairwise alignments between a single sequence and each sequence withing a list of sequences
    
    Object attributes:
        - target_record (str): the targeted sequence, each alignment would be done with respect to it
        - breeds_records (list[str, ...]): a list of sequences, each sequence would be aligned with the targeted sequence seperetaly. 
        - scores (list[float, ...]): the scores of pairwise alignments


    Object methods:
        - __init__(self, graph): Create a DistributionCalc object
        - fit_gumbel_distribution(self): Computes the Lambda.
        - compute_e_value(self, score, lambda_, K): Computes the E value.
        - compute_p_value(self, e_value): Computes a single p_value.
        - compute_p_values(self): Computes the p_value from provided scores list.
        - p_values_get(self): Returns the calculated p value.
    """

    def __init__(self, target_record, breeds_record, scores):

        self.target_record = target_record
        self.breeds_record = breeds_record

        self.scores = scores
        self.p_values = np.zeros(len(scores))

        self.compute_p_values()


    def fit_gumbel_distribution(self):
        """ Computes the Lambda. """
        shape, loc, scale = genextreme.fit(-self.scores) 
        lambda_const = 1 / scale 
        K = 0.1  
        return lambda_const, K
    

    def compute_e_value(self, score, lambda_, K):
        """ Computes the E value. """
        M, N = len(self.target_record.seq), len(self.breeds_record[0].seq)  # Sequence lengths
        return K * M * N * np.exp(-lambda_ * score)
    

    def compute_p_value(self, e_value):
        """ Computes a single p_value. """
        return 1 - np.exp(-e_value)
    

    def compute_p_values(self):
        """ Computes the p_value from provided scores list. """
        lambda_const, K = self.fit_gumbel_distribution()
        for i, score in enumerate(self.scores):
            e_value = self.compute_e_value(score, lambda_const, K)
            self.p_values[i] = self.compute_p_value(e_value)


    def p_values_get(self):
        """ Returns the calculated p value. """
        return self.p_values