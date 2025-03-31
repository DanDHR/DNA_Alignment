_This code identifies the most similar sequence between a mystery dna sequence of a dog breed and a list of dna sequences from another 99 dog breeds. Futhermore, the code will provide the user with additional information, in the form of:_<br />
_1) Shell output that details the statistics of the sequence similarity scores._<br />
_2) A plot of the phylogenetic tree of those dog breeds._<br /> 

# DNA_Alignment

# To run the code: <br /> 
1. Clone repo from github: `git clone https://github.com/DanDHR/DNA_Alignment` <br />
2. To ensure no dependence clashes, use: <br />
        `python -m venv test_env` - Create a new python environment<br />
        `source test_env/bin/activate` - Activate the new environment<br />
       `pip install -r requirements.txt` - Install the all dependences<br />
3. `cd` to `DNA_Alignment` if not already there and run  `python3 Main.py`. <br />




# Details on phylogeny
For the phylogenetic tree, a multi sequence alignment is needed. For that purpose the use of Clustal Omega (https://www.ebi.ac.uk/jdispatcher/msa/clustalo) was required. 
The input parameters are shown bellow:
![image](https://github.com/user-attachments/assets/2177b76c-5d20-4dfe-ae60-5640d606db5b)


















