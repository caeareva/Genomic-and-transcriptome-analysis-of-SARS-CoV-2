#!/usr/bin/env python3

##################################################################################################
#	
#	Author: Carlos Arevalo (caeareva)
#
#	File: seqTranslator.py 
#	Executable: python seqTranslator.py < sars2FilProt.fa > sars2ProteinSeq.fa 
#				python seqTranslator.py <  btRsFiltProt.fa > btRsProteinSeq.fa
#				python seqTranslator.py < sars4231FilProt.fa > sars4231ProteinSeq.fa
#				python seqTranslator.py < sarsG13FilProt.fa > sarsG13ProteinSeq.fa
#				python seqTranslator.py < cvUrbaniFilProt.fa > coronavirusUProteinSeq.fa
#
#	Required module: FastAreader
#	Purpose: translate multiple RNA and DNA sequences to single letter amino acid sequences
#	Condition(s): None
#
##################################################################################################

'''
In this class, we define objects to read FastaA files inside the sequenceAnalysis program.
'''

import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
           
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:

            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

"""
Class translated RNA and DNA to protein sequence
"""
class translateSeq:
	"""This calss translate DNA and RNA sequences to amino acid sequences"""

	def validateDNA(dnaSequence):
		"""
		Validates DNA sequence. Checks for correct nucleotides and returns 
		True if sequence is valid.
		"""
		seq = dnaSequence.upper()
		# replace all unwated characters in sequence
		seq = seq.replace("_", "").replace(" ", "").replace("-", "").replace("None", "")
		validSeq = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
		# if sequence is valid return it
		if validSeq == len(seq):
			return True
		else:
			return False 

	def validateAmino(aaSequence):
		"""
		Validates amino acid sequence. Checks for correct amino acids and returns True
		if sequence is valid.
		"""
		seq = aaSequence
		# remove unwated character from protein sequence
		seq = seq.replace("_", "").replace(" ", "").replace("-", "")
		validAminos =	seq.count("A") + seq.count("V") + seq.count("I") + seq.count("L") + \
						seq.count("M") + seq.count("F") + seq.count("Y") + seq.count("W") + \
						seq.count("S") + seq.count("T") + seq.count("N") + seq.count("Q") + \
						seq.count("C") + seq.count("G") + seq.count("P") + seq.count("R") + \
						seq.count("H") + seq.count("K") + seq.count("D") + seq.count("E") + \
						seq.count("_") 
		# if protein sequence is valid, return it
		if validAminos == len(seq):
			return True
		else:
			return False 

	def translateCodon(codon):
		"""Translate a codon into a single letter amino acid"""

		rnaCodonTable = {
        	# Second Base
        	# U             C             A             G
        	# U
        	'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',
        	'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
        	'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',
        	'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',
        	# C
        	'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',
        	'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        	'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        	'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
        	# A
        	'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',
        	'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        	'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        	'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
        	# G
        	'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',
        	'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        	'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        	'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
    		}

    	# converts RNA codon table to DNA codon table
		dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

		# stores dictionaries 
		context = dict(rnaCodonTable.items() | dnaCodonTable.items()) 
   		
   		# if codon is found in dictionary return it
		if codon in context:
			return context[codon]
		else:
			return None

	def translator(dnaSequence, start = 0):
		"""Translate a DNA sequence into amino acid sequence (protein).
		Return a putative protein sequence"""
		
		fastaFile = FastAreader()
		dnaSequence = fastaFile.readFasta()
		for header, sequence in fastaFile.readFasta():
			cleanSeq = sequence.replace("None", "")
			# validate incoming sequence
			assert translateSeq.validateDNA(cleanSeq), """Invalid DNA sequence"""
			print(">{}".format(header))
			aaSeq = " " 
			# translate valid sequences only
			for position in range(start, len(cleanSeq)-2, 3):
				codon = cleanSeq[position: position + 3]
				protein = translateSeq.translateCodon(codon)
				aaSeq += protein
			
			print(aaSeq.replace("-", "") + "\n")
			
def main():
	"""Call functions and return putative proteins from ORFs"""
	aminoSeq = translateSeq()
	print(aminoSeq.translator())

if __name__ == "__main__":
    main()





