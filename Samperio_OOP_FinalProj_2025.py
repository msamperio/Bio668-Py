## RENAME this file YourLastName_OOP_FinalProject_2023.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes overlapping kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT

#sequence.reverse()

### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames
    #Original sequence shifted by 0, 1, and 2 bases
    
### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object. 


import re
import doctest

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}


class Seq:

    def __init__(self, sequence='', gene= '', species='', kmers=[]):
        """
        Initializes the sequence, gene, species, and kmers list. Ensures sequence 
        is upper case and stripped.
        
        >>> s = Seq("ATGCCGTA", "Gene1", "Homo sapiens")
        >>> s.sequence
        'ATGCCGTA'
        >>> s.gene
        'Gene1'
        >>> s.species
        'Homo sapiens'
        """
        
        self.sequence=sequence.strip().upper()
        self.gene=gene
        self.species=species
        self.kmers=[] #needs to be empty list
        
    def __str__(self):
        """
       Returns the gene and species information in a string.
       
       >>> s = Seq("ATGC", "geneA", "Homo sapiens")
       >>> s.print_record()
       'Homo sapiens geneA: ATGC'
       """
        return self.sequence

    def print_record(self):
        """
        Returns the gene and species information in a string.
        
        >>> s = Seq("ATGC", "geneA", "Homo sapiens")
        >>> s.print_record()
        'Homo sapiens geneA: ATGC'
        """
        return self.species + " " + self.gene + ": " + self.sequence

    def make_kmers(self, k=3):
        """
       Generates k-mers from the sequence. Default parameter is 3.
       
       >>> s = Seq("ATGC", "geneA", "Homo sapiens")
       >>> s.make_kmers(2)
       ['AT', 'TG', 'GC']
       """
        for i in range(0,len(self.sequence)):
            kmer=self.sequence[i:i+k] #in case kmer parameter is changed.
            #print(kmer)
            #make it so kmers length are only equal to k 
            if len(kmer)==k:
                self.kmers.append(kmer)
        return self.kmers
  
    def fasta(self):
        """
        Return the sequence in FASTA format.
        
        >>> s = Seq("ATGCCGTA", "Gene1", "Homo sapiens")
        >>> print(s.fasta())
        >Homo sapiens Gene1
        ATGCCGTA
        """
        return ">"+ self.species+" "+ self.gene  + "\n"+self.sequence
    
class DNA(Seq):
    
    def __init__(self, sequence='', gene= '', species='', geneid='',kmers=[]):
        """
        Initialize the DNA sequence, including geneid, and introduces N to sequence.
        
        >>> d = DNA("GAxTCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> d.sequence
        'GANTCTC'
        >>> d.geneid
        'GeneID123'
        """
        super().__init__(sequence,gene,species,kmers)
        self.sequence=re.sub('[^ATGCU]','N',self.sequence)
        self.geneid=geneid
        
    def __str__(self):
        """
        Return the string representation of the DNA sequence.
        
        >>> d = DNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> str(d)
        'GATCTC'
        """
        return super().__str__()
    
    def analysis(self):
        """
        Calculate GC content of the DNA sequence.
        
        >>> d = DNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> d.analysis()
        3
        """
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
        """
        Print DNA record including geneid, species, and sequence.
        
        >>> d = DNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> d.print_info()
        'GeneID123 Homo sapiens Gene1: GATCTC'
        """
        return self.geneid +" "+ self.species + " " + self.gene + ": " + self.sequence  

    def reverse_complement(self):
        """
        Get the reverse complement of the DNA sequence.
        
        >>> d = DNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> d.reverse_complement()
        'GAGATC'
        """
        rseq = self.sequence[::-1]
        d = {"A": "T", "T": "A", "C": "G", "G": "C","N":"N"}
        compseq = ""
        for base in rseq:
            compseq += d[base]
        return compseq
        
    def six_frames(self):
        """
        Generate the six reading frames (3 forward, 3 reverse).
        
        >>> d = DNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> d.six_frames()
        ['GATCTC', 'ATCTC', 'TCTC', 'GAGATC', 'AGATC', 'GATC']
        """
        frames=[]
        for i in range(3): #will stop at third character in string
            frames.append(self.sequence[i:])#need to slice or else won't work
        
        rev_comp = self.reverse_complement()
        for i in range(3):  
            frames.append(rev_comp[i:])
        return frames


class RNA(DNA):

    def __init__(self, sequence='', gene= '', species='', geneid='',kmers=[],codons=[]):
        """
        Initialize the RNA sequence by converting DNA sequence to RNA.
        
        >>> r = RNA("GATCTC", "Gene1", "Homo sapiens", "GeneID123")
        >>> r.sequence
        'GAUCUC'
        """
        super().__init__(sequence,gene,species,geneid,kmers)
        self.codons=[]
        self.sequence=re.sub('T','U',self.sequence)
        
    def make_codons(self): 
       """
       Create codons from the RNA sequence and store them in self.codons.
       >>> r = RNA("GAUCUCGAA", "Gene1", "Homo sapiens", "GeneID123")
       >>> r.make_codons()
       ['GAU', 'CUC', 'GAA']
       >>> r.codons
       ['GAU', 'CUC', 'GAA']
       """
       self.codons = []
       for i in range(0, len(self.sequence), 3):
           codon = self.sequence[i:i+3]  
           if len(codon) == 3:  
               self.codons.append(codon) 
       return self.codons
   
    def translate(self):
        """
        Translate the RNA codons into a protein sequence.
    
        >>> r = RNA("GAUCUC", "Gene1", "Homo sapiens", "GeneID123")
        >>> r.make_codons()
        ['GAU', 'CUC']
        >>> r.translate()
        'DL'
        """
        proteinlist=[] #output should be string
        for codon in self.codons:
            if codon in standard_code:  
                proteinlist.append(standard_code[codon])#append only works with lists 
            else:
                pass
        proteinfinal = "".join(proteinlist) #"joins" contents of a list into a string
        return proteinfinal
            
class Protein(Seq):

    def __init__(self, sequence='', gene= '', species='', kmers=[]):
        """
        Initialize the protein sequence and change non letter characters to X.
    
        >>> p = Protein("AC1DE", "Gene1", "Homo sapiens")
        >>> p.sequence
        'ACXDE'
        """
        super().__init__(sequence,gene,species,kmers)
        self.sequence = re.sub('[^A-Z]', 'X', self.sequence)#caret ^ match anything that is NOT a letter
        
    def total_hydro(self):
        """
        Calculate total hydrophobicity based on Kyte-Doolittle scale.
    
        >>> testp=Protein('VIKING','test','unknown',999)
        >>> x=testp.total_hydro()
        >>> print(x)
        5.4
        """
        hvalues = []  
        
        for aa in self.sequence:  
            if aa in kyte_doolittle:  
                hvalues.append(kyte_doolittle[aa])  
        
        totalh = sum(hvalues)  
        return totalh 

    def mol_weight(self):
        """
        Calculate the molecular weight of the protein.

        >>> testp=Protein('VIKING','test','unknown',999)
        >>> m=testp.mol_weight()
        >>> print(m)
        732.87
        """
        wvalues = []  
        
        for aa in self.sequence:  
            if aa in aa_mol_weights:  
                wvalues.append(aa_mol_weights[aa])  
        
        totalw = sum(wvalues)  
        return totalw


    

x=DNA("G","tmp","m",000)



if __name__ == "__main__":
    doctest.testmod(verbose=True)





