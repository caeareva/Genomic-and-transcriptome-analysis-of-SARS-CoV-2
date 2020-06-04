# !/usr/bin/env python3

#################################################################################################
#
#   Author: Carlos Arevalo (caeareva)
#  
#   File: getCodingSeq.py
#   Executable: python getCodingSeq.py < sars2FilteredSeq.fa > sars2PutativeSeq.fa
#               python getCodingSeq.py < coronavirusBtRsSeq.fa > btRsPutativeSeq.fa
#               python getCodingSeq.py < sarsG13Seq.fa > sarsG13PutativeSeq.fa
#               python getCodingSeq.py < sarsRs4231Seq.fa > sarsRs4231PutativeSeq.fa
#               python getCodingSeq.py < coronavirusUrbaniSeq.fa > cvUrbaniPutativeSeq.fa
#
#   Required module: FastAreader
#   Purpose: Obtain putative proteins (codin sequences) from each gene or fasta file
#   Condition(s): None 
#
##################################################################################################

import sys

'''
In this class, we define objects to read FastaA files inside the sequenceAnalysis program.
'''
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

        yield header, sequence

'''
Program get coding sequences from gene fasta files.
In Eukaryotes, a coding sequence is defined from the first start codon to the first stop codon in an ORF.
'''

def main():
    '''
    Function read the fasta file and defines a putative sequence. Program defines putative sequences 
    from the fisrt start codon to the first stop codon. Program outputs a fasta file containing putative
    proteins only.
    '''
    myReader = FastAreader()
    for header, sequence in myReader.readFasta():
        print('>{}'.format(header))
        sequence = sequence.upper()
        sequence = sequence.replace('_', '')
        sequence = sequence.replace('\n', '')
        sequence = sequence.replace('\r', '')
        for start in range(0, len(sequence), 3):
            codon = sequence[start: start + 3]
            if codon == 'ATG':
                startPosition = sequence.find('ATG')
                proteinSeq = sequence[startPosition:]
                
        print(proteinSeq + '\n')

if __name__ == "__main__":
    main()
                        
