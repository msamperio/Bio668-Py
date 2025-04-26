#------------------------------------------------------------------------------------------------
#PYTHON FINAL PRACTICAL
#------------------------------------------------------------------------------------------------
#READ CAREFULLY AND THOROUGHLY. TEST THE FUNCTIONS THOROUGHLY.
#WHEN YOU ARE DONE,
#(1) SAVE THE FILE NAME WITH YOUR LAST NAME PYTHON FILE Ex: KelleyPythonExam2_2024.py
#(2) MAKE SURE THE FILE RUNS - NO ERRORS - EVEN IF YOU DIDN'T GET EVERYTHING CORRECT!
#(3) WHEN I RUN THE FILE IN PYTHON, NOTHING SHOULD BE PRINTED OUT (POINTS WILL BE SUBTRACTED)
#(4) SUBMIT YOUR EXAM

#------------------------------------------------------------------------------------------------
#== FUNCTION 1 ==
#Write a function called blosum that takes 2 parameters, observed and expected
# values and returns the blosum score. 
# The equation: blosum_score=2*log2(observed/expected). 
import math

def blosum(obs,exp):
    math1=math.log2(obs/exp)
    blosum_score=2*math1
    return blosum_score

#Test with obs and exp
obs=0.8
exp=0.1
#print(blosum(obs, exp))
#Should print: 6.0

#------------------------------------------------------------------------------------------------
#== FUNCTION 2 ==
#Write a function called check_start_codon that takes 1 parameter,
# a list of dna or rna sequences and returns a list of all the sequences
# in which the first 3 letters are a start codon: "ATG" or "AUG"
"""
import re
def check_start_codon(seq):
    match=[""]
    match1=re.match("AUG",seq)
    match2=re.match("ATG",seq)
    match.append(match1)
    return match1, match2
"""
"""
def check_start_codon(seq):
    list=[""]
    for b in range(0,len(seq),3):
        codon=seq[b:b+3]
        if codon== "AUG" or codon== "ATG":
            list.append(codon)
        return list
"""      
#Test check_start_codon with the following:
mixed_seqs=["AUGGAGACA", "ATGTGC", "ATCGGAT","AUGGGACAG","ATATAT"]
#this=check_start_codon(mixed_seqs)
#print(this)
#Should return: ['AUGGAGACA', 'ATGTGC', 'AUGGGACAG']

#------------------------------------------------------------------------------------------------
#== FUNCTION 3 ==
#Write a function called alpha_total that takes 1 parameter,a protein string
# and uses the cf dictionary which has P(a) scores to return the total P(a) score
# for the entire sequence. 
# NOTES: (1) The function is partially completed
#        (2) It should work with clean and 'dirty' sequence (see below)

cf={"R":98,"K":114,"N":67,"D":101,"Q":111,"E":151,"H":100,
        "P":57,"Y":69,"W":108,"S":77,"T":83,"G":57,"A":142,"M":145,
        "C":70,"F":113,"L":121,"V":106,"I":108}
"""
def alpha_total(prot):
    count=0
    clean=prot.strip()
    cleaner=clean.upper()
    for aa in cleaner:
        pa=cf[aa]
        total=
    return total
         
test_prot1 = "QYVED"
#Test the function with test_prot1
print(alpha_total(test_prot1))
test_prot2 = "\t QyVEd\n"
#Test the function with test_prot2
print(alpha_total(test_prot2))

"""