import random

def countStarts( DNA ):
    ''' Given a string of A, T, C, and G's, returns the number of
        codons that are ATG.'''
    counter = 0
    for i in range(len(DNA)):
        if DNA[i:i+3] == 'ATG':
            counter = counter + 1
    return counter

def genString( GCcontent, length ):
    ''' Given the GC content (between 0 and 1) and the desired length,
        returns a randomly generated DNA string with the given GC
        content and length.'''
    output = ''
    for i in range(length):
        myRandomNumber = random.random()
        if myRandomNumber < GCcontent:
            nuc = random.choice(['G', 'C'])
        else:
            nuc = random.choice(['A', 'T'])
        output = output + nuc
    return output


def startCodons( GCcontent, length, trials ):
    ''' Given the desired GC content (between 0 and 1), the length of
      the random DNA strings, and the number of trials to perform,
      returns the average fraction of codons that are start codons.'''
    total = 0
    for r in range(trials):
        DNAseq = genString(GCcontent, length)
        total = total + countStarts(DNAseq)
    return 1.0 * total / (trials * (length-2))