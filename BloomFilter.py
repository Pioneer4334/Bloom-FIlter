import math
import random
from bitarray import bitarray

class BloomFilter:
    """
    constructor that takes as input a list of n DNA strings (of the virus) and a desirable false positive error probability p
    """
    def __init__(self, viralGenome, expFalsePositive):
        self.p = expFalsePositive
        self.viralGenome = viralGenome
        self.n = len(viralGenome)
        self.l = len(viralGenome[0])
        self.m = int(-self.n*math.log(expFalsePositive)/((math.log(2))**2))
        self.k = int((self.m/self.n)*math.log(2))
        self.bf = bitarray('0'*self.m)   
    
    """
    generates k random hash functions for DNA strings of length L
    """
    def fn_Hash(self, dna):
        dictDNA = { "A":0, "C":1, "G":2, "T":3 }
        l=len(dna)
        hash_list = list()
        if(l == self.l):
            for i in range(self.k): 
                h = 2166136261
                random.seed(i+1)
                for j in range(l):
                    h = (random.randint(1, self.k) + h * 16777619) ^ (dictDNA[(dna[j].upper())] + 97)
                hash_list.append((h) % self.m)
        return hash_list
    
    """
    inserts a DNA string into the Bloom filter
    """
    def fn_Insert(self, dna):
        hashedIndices = self.fn_Hash(dna)
        for hashedIndex in hashedIndices:
            self.bf[hashedIndex] = True
    
    """
    checks if an unknown DNA string is a viral DNA
    """
    def fn_CheckDNA(self, dna):
        hashedIndices = self.fn_Hash(dna)
        for hashedIndex in hashedIndices:
            if(self.bf[hashedIndex] == False):
                return False
        return True

"""
generates a random DNA of length L.
"""
def fn_GenerateDNA(length):
    dna=""
    random.seed(None)
    for i in range(length):
        dna +=random.choice(['A', 'C', 'G', 'T'])
    return dna

"""
evaluation function that generates a viral genome containing number of random DNA strings of specified length and test to see if the bloom filter works
"""   
def fn_evaluate(noOfDNA, lengthOfDNA, noOfTestDNA, expectedProb):
    viralDNA = list()
    for i in range(noOfDNA):
        viralDNA.append(fn_GenerateDNA(lengthOfDNA))

    bloom = BloomFilter(viralDNA, expectedProb)
    for i in range(bloom.n):
        bloom.fn_Insert(bloom.viralGenome[i])

    fp = 0
    for i in range(noOfTestDNA):
        unknownDNA = fn_GenerateDNA(lengthOfDNA)
        flagPresence = bloom.fn_CheckDNA(unknownDNA)
        print("\033[93m" if flagPresence else "\033[92m", "The DNA", unknownDNA, "is", "a virus" if flagPresence else "not a virus", "\033[0m")
        if(flagPresence and unknownDNA not in viralDNA):
            fp += 1
            print("\033[91m !!!False Positive!!! \033[0m")
        elif (not flagPresence and unknownDNA not in viralDNA):
            print("\033[94m !!!True Negative!!! \033[0m")
        elif (flagPresence and unknownDNA in viralDNA):
            print("\033[94m !!!True Positive!!! \033[0m")
        else:
            print("\033[94m !!!False Negative!!! \033[0m")

    return fp/noOfTestDNA

expectedFP=0.01
print("Expected False Positive Probability:", expectedFP, "Obtained False Positive Probability", fn_evaluate(noOfDNA=10000, lengthOfDNA=100, noOfTestDNA=1000, expectedProb=expectedFP))
    
    