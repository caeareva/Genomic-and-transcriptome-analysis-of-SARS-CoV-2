#!/usr/bin/env python3

#################################################################################################
#
#   Author: Carlos Arevalo (caeareva)
#
#   File: filterORFs.py
#   Executable: python filterORFs.py < sars2PutativeSeq.fa > sars2FilProt.fa
#			    python filterORFs.py < sarsG13PutativeSeq.fa > sarsG13FilProt.fa
#               python filterORFs.py < sarsRs4231PutativeSeq.fa > sars4231FilProt.fa
#               python filterORFs.py < btRsPutativeSeq.fa > btRsFiltProt.fa
#               python filterORFs.py < cvUrbaniPutativeSeq.fa > cvUrbaniFilProt.fa
#
#   Required module: FastAreader
#   Purpose: match putative protein with scientific name and ignored nonfunctional proteins
#   Condition(s): Dictinary need to be expanded due that some proteins need to be experimentally
#                 study to determine their functionality (focused on proteins already validated).
#
##################################################################################################

'''
In this class, we define objects to read FastaA files inside the sequenceAnalysis program.
'''
import sys
class FastAreader:
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


'''
In this class, we defined ORFs coordinates and find their ORFs names in dictionaries.
'''
class nameFinder:
    '''
    Read dictionaries and match putative proteins with scientific names. Nonfunctional 
    proteins are ignored. Output a file containing proteins of interest.
    '''

    covidORFs = {'SARSCoV2 ORF1a': 'NC_045512.2:266-13483', 'SARSCoV2 ORF1ab': 'NC_045512.2:266-21555', 
                'SARSCoV2 S': 'NC_045512.2:21562-25384', 'SARSCoV2 S': 'NC_045512.2:21536-25384',
                'SARSCoV2 ORF3a': 'NC_045512.2:25393-26220', 'SARSCoV2 E': 'NC_045512.2:26245-26472', 
                'SARSCoV2 M': 'NC_045512.2:26523-27191', 'SARSCoV2 ORF6': 'NC_045512.2:27202-27387', 
                'SARSCoV2 ORF7a': 'NC_045512.2:27394-27759', 'SARSCoV2 ORF7b' : 'NC_045512.2:27756-27887',
                'SARSCoV2 ORF8': 'NC_045512.2:27894-28259', 'SARSCoV2 N': 'NC_045512.2:28275-29533', 
                'SARSCoV2 ORF10': 'NC_045512.2:29558-29674'}

    ratG13 = {'ratG13 ORF1ab': 'MN996532.1:251-21537', 'ratG13 S': 'MN996532.1:21545-25354',
              'ratG13 NS3': 'MN996532.1:25363-26190', 'ratG13 N': 'MN996532.1:28240-29499',
              'ratG13 NS8': 'MN996532.1:27860-28225', 'ratG13 M': 'MN996532.1:26493-27158',
              'ratG13 NS7a': 'MN996532.1:27360-27725', 'ratG13 E': 'MN996532.1:26215-26442',
              'ratG13 NS6': 'MN996532.1:27169-27354', 'ratG13 NS7b': 'MN996532.1:27722-27853',
              'ratG13 ORF1a': 'MN996532.1:251-13465' }

    coronavirusUrbani = {'CoV_U S2': 'AY278741.1:648-1252', 'CoV_U NS8': 'AY278741.1:7-79',
                         'CoV_U M': 'AY278741.1:3-220', 'CoV_U ORF3b': 'AY278741.1:1-153',
                         'CoV_U S': 'AY278741.1:21492-25259', 'CoV_U S receptor binding': 'AY278741.1:318-569',
                         'CoV_U M': 'AY278741.1:26398-27063','CoV_U PP1a': 'AY278741.1:265-13413' }


    rs4231 = {'Rs4231 S2': 'KY417146.1:648-1252', 'Rs4231 NS8': 'KY417146.1:1-118',
              'Rs4231 S receptor binding': 'KY417146.1:318-569', 'Rs4231 S': 'KY417146.1:21493-25260',
              'Rs4231 M': 'KY417146.1:3-220', 'Rs4231 ORF3b': 'KY417146.1:1-114',
              'Rs4231 ORF6': 'KY417146.1:27075-27266', 'Rs4231 ORF7b': 'KY417146.1:27639-27773',
              'Rs4231 ORF8': 'KY417146.1:27780-28145', 'Rs4231 ORF7a': 'KY417146.1:27274-27642',
              'Rs4231 ORF3a': 'KY417146.1:25269-26093', 'Rs4231 ORF1a': 'KY417146.1:265-13413'}

    BtRsBetaCoV = {'BtRsBetaCoV ORF1ab': 'MK211376.1:264-21484', 'BtRsBetaCoV spike glycoprotein': 'MK211376.1:21491-25261',
                   'BtRsBetaCoV Matrix': 'MK211376.1:26400-27065', 'BtRsBetaCoV S': 'MK211376.1:21491-25261', 
                   'BtRsBetaCoV nucleocapsid': 'MK211376.1:28683-29951'}
        
    # provides access to both dictionaries in single object 
    context = dict(covidORFs.items() | ratG13.items() | rs4231.items() | BtRsBetaCoV.items() | coronavirusUrbani.items())

    def findName(self):
        '''
        Read thoughout dictionaries and matches genes
        '''
        orfDict = self.context.items() # access dictionaries items
        inFasta = FastAreader()
        for fasta in inFasta.readFasta():
            header, sequence = fasta
            pass
            
            # acces genes fasta files
            for key, value in orfDict: 
                if header == value:
                    # matches gene with scientific name
                    header = key
                    print(">{}".format(key)) 
                    print(sequence + '\n')
                    pass
                else:
                   pass
             
def main():
    '''
    Call functions and output file
    '''
    myName = nameFinder()
    output = myName.findName()
    print(output) 

if __name__ == "__main__":
    main()


