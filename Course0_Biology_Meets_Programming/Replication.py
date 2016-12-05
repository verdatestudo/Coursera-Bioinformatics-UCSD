'''
Replication.py
Biology Meets Programming: Bioinformatics for Beginners - UCSD - Week 1

Last Updated: 2016-Dec-03
First Created: 2016-Dec-03
Python 3.5
Chris
'''

def PatternCount(Pattern, Text):
    '''
    Takes a string 'Pattern' and returns an integer of the number of occurances in a 'Text' string.
    e.g 'AB', 'ABABCCCAB' returns 3.
    '''
    return len([x for x in range(len(Text) - len(Pattern) + 1) if Text[x:x+len(Pattern)] == Pattern])

def CountDict(Text, k):
    '''
    Takes a string 'Text' and integer 'k'.
    Returns a dict of how many times each k-mer of a sliding window occurs.
    e.g 'Text' = 'CGATATA'. dict[0] = 1 (#CGA), dict[1] = 1 (#GAT), dict[2] = 2 (#ATA) and so on.
    '''
    return {i: PatternCount(Text[i:i+k], Text) for i in range(len(Text) - k+1)}

def FrequentWords(Text, k):
    '''
    'Text' string and k integer. Returns a set of strings.
    Takes a text string and finds the most frequent patterns of k length.
    '''
    # have to calc max separately in case there are multiple keys of the same max value.
    Count = CountDict(Text, k)
    m = max(Count.values())

    return set(Text[i:i+k] for i in Count if Count[i] == m)

def ReverseComplement(Pattern):
    '''
    Reverse Complement Problem: Find the reverse complement of a DNA string.
    Input: A DNA string Pattern.
    Output: The reverse complement of Pattern.
    '''
    # replace each letter with it's complement and reverse the string (due to 5 prime and 3 prime DNA)
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C':'G'}
    return ''.join([comp[x] for x in Pattern])[::-1]

def PatternMatching(Pattern, Genome):
    '''
    # Input:  Two strings, Pattern and Genome.
    # Output: A list containing all starting positions where Pattern appears as a substring of Genome.
    '''
    return [x for x in range(len(Genome) - len(Pattern) + 1) if Genome[x:x+len(Pattern)] == Pattern]

def questions():
    '''
    Questions.
    '''

    print(CountDict('CCGAACACCCGTACACCGAACACCACACCACACCTTGCACACCACACCTACACCACACACCACACCGGACACCCACACCCACACCACGAACACCGAGAGTACACCTA', 5))

    vibrio_cholerae_1 = 'ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTAT'\
    'CTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTC'\
    'CAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTT'\
    'GACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACAT'\
    'GCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTG'\
    'ATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'

    Pattern = 'TGATCA'
    print(PatternCount(Pattern, vibrio_cholerae))

    words = FrequentWords("GATCCAGATCCCCATAC", 2)
    print(words)

    print(FrequentWords('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4))

    print(FrequentWords(vibrio_cholerae_1, 10))

    print(ReverseComplement("ATGATCAAG"))
    print(ReverseComplement("CTTGATCAT"))
    print(ReverseComplement("TCTTGATCA"))
    print(ReverseComplement("CTCTTGATC"))

    print(PatternMatching('ATAT', 'GATATATGCATATACTT'))

    thermotoga_petrophila = 'AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTT'\
    'GTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTT'\
    'AGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCG'\
    'TTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAA'\
    'CCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA'

    count1 = PatternCount("ATGATCAAG", thermotoga_petrophila)
    count2 = PatternCount("CTTGATCAT", thermotoga_petrophila)
    print(count1+count2) # equals 0, Thus, different bacteria may use different DnaA boxes as “hidden messages” to the DnaA protein.

    print(PatternCount("CGCG", "CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC"))
    print(FrequentWords("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3))
    print(ReverseComplement("TTGTGTC"))
    print(True or not False and False)
    
#questions()
