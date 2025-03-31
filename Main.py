from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Pairwise_Alignment import PairwiseAlignment
from DistributionCalc import DistributionCalc
from Phylogeny import PhylogeneticTreeBuilder
import matplotlib.pyplot as plt


def fasta_file_reader(file_path) -> list:
    """
    Reads a fasta file and returns all sequences within it as a list
    
    Args:
        - file_path (str): the path and the file that is going to be read

    Return:
        - list: list of the sequences within the read file
    """
    #Read the file and set its type as a list
    record = list(SeqIO.parse(file_path, "fasta"))
    return record   


def find_best_pairwise_alignment(target_record, breeds_record, difference = "") -> list:
    """
    Finds the sequence from all of the breeds records that is closest to the targeted record based on their alignment scores. Furthermore,
    it calculates their difference as a percent.

    Args:
        - target_record (str): the record that we are interested in
        - breeds_record (list[str,...]): a list of sequences from dog breeds
        - difference (str): if set to "calculate" it will calculate the percentage of difference between the best scored sequence from the
            breeds_record and the targeted record. 
            !!!The default is set to not calculate this, since it takes too much recourses and it takes a long time to process that!!!
    Return:
        - list: a list of all the alignment scores between each of the breeds records and the targeted record
    """
    #Initialise PairwiseAlignment class. It will save an instance of each of those records
    pa = PairwiseAlignment(target_record, breeds_record)
    #Calculate all alignments between the targeted record and the breeds records
    pa.all_pairwise_alignments()
    #Get the location of the best aligned sequence within the breeds records
    best_index = pa.best_index_get()

    #Print results
    print("\nClosest match:", breeds_records[best_index].id)
    print("Sequence:", breeds_records[best_index].seq)

    #Calculate the difference as a percentage between the targeted record and the best aligned sequence to it from the breeds records
    if difference == "calculate":
        #Since to calculate that, we need the actual alignments, it takes too long to process it
        alignment = pa.single_pairwise_alignment(3, best_index)
        print(f"Difference percentage: {pa.difference_percent(alignment):.2f}%")

    return pa.scores_get()      #Return the scores of all alignments. It is needed later.


def records_distribution(target_record, breeds_record, scores) ->list:
    """
    Calculate the distribution of sequence scores from sequence alignments

    Args:
        - target_record (str): the record that we are interested in
        - breeds_record (list[str,...]): a list of sequences from dog breeds
        - scores (list[float, ...]): the scores of the pairwise alignments between the target record and each of the breeds records

    Return:
        -list: for each score, calculate its p-value    
    """
    distr = DistributionCalc(target_record, breeds_record, scores)          #Initialise the DistributionCalc to calculate the p-values from the scores
    return distr.p_values_get()                                             #Return the p-values


if __name__ == '__main__':
    #Load the target and the list of dog breeds dna sequences
    breeds_records = fasta_file_reader("project_dog_dna/dog_breeds.fa")
    target_record = fasta_file_reader("project_dog_dna/mystery.fa")
    #Find the best pairwise alignment and get their scores
    scores = find_best_pairwise_alignment(target_record[0],breeds_records)          #Main goal


    #Calculate the distribution of sequence scores from sequence alignments
    print(records_distribution(target_record[0], breeds_records, scores))           #Stretch goal 1


    #Load a Multiple Sequence Alignment
    alignment = AlignIO.read("external_data/dog_breeds_alignment.fa", "fasta")
    #Create a phylogenetic tree
    phylo_tree = PhylogeneticTreeBuilder(alignment)                                 #Stretch goal 2
