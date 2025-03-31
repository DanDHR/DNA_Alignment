import numpy as np
import unittest
from scipy.stats import genextreme 


##The code bellow is inspired by the statistics of sequence similarity scores done by Blast: https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html


class DistributionCalc(unittest.TestCase):
    """
    A class with the purpose of calculating the distribution of sequence scores from sequence alignments
    
    Object attributes:
        - target_record (str): the targeted sequence, each alignment would be done with respect to it
        - breeds_records (list[str, ...]): a list of sequences, each sequence would be aligned with the targeted sequence seperetaly. 
        - scores (list[float, ...]): the scores of pairwise alignments


    Object methods:
        - __init__(self, graph): Create a DistributionCalc object
        - fit_gumbel_distribution(self): Initialises both the lambda and K constant values.
        - compute_e_value(self, score, lambda_, K): Computes the E value.
        - compute_p_value(self, e_value): Computes a single p_value.
        - compute_p_values(self): Computes the p_value from provided scores list.
        - p_values_get(self): Returns the calculated p value.
    """

    def __init__(self, target_record, breeds_record, scores):
        #We need both records only for their lengths
        self.target_record = target_record
        self.breeds_record = breeds_record
        #Get the scores of all alignments
        self.scores = scores
        #Initialise a place holder array for p-values
        self.p_values = np.zeros(len(scores))
        #Compute all the p-values
        self.compute_p_values()


    def fit_gumbel_distribution(self) -> tuple[float, float]:
        """ 
        Initialises both the lambda and K constant values. 
        """
        shape, loc, scale = genextreme.fit(-self.scores)    

        #Compute lambda
        lambda_const = 1 / scale 

        #Set K
        K = 0.1                 
        return lambda_const, K
    

    def compute_e_value(self, score, lambda_, K) -> float:
        """ 
        Computes the E value.

        Args:
            - score (float): the provided score for an alignment
            - lambda_ (float): lambda constant 
            - K (float): K constant
            
        Return:
            - float: the computed E value
        
        """
        M, N = len(self.target_record.seq), len(self.breeds_record[0].seq)  #Sequence lengths
        return K * M * N * np.exp(-lambda_ * score)                         #Return the E value
    

    def compute_p_value(self, e_value) -> float:
        """ 
        Computes a single p_value.

        Args:
            - e_value (float): the previously computed E value

        Return:
            - float: the coputed p-value
        
        """
        return 1 - np.exp(-e_value)
    

    def compute_p_values(self) -> None:
        """ 
        Computes the p_values from provided scores list. 
        """

        lambda_const, K = self.fit_gumbel_distribution()                            #Initialise lambda and K constants

        #iterate through all scores
        for i, score in enumerate(self.scores):
            e_value = self.compute_e_value(score, lambda_const, K)                  #Calculate the E value
            self.p_values[i] = self.compute_p_value(e_value)                        #Calculate the p-value


    def p_values_get(self) -> float:
        """ 
        Returns the calculated p value. 
        """
        return self.p_values


    def p_values_get(self):
        """ Returns the calculated p value. """
        return self.p_values
