#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 10:03:04 2020

@author: pi
"""

'''
Replication.py contains functions to find patterns in DNA sequences.

From a sequence, these functions can:
    - Search for a kmer.
    - Find frequent kmers.
    - Generate the reverse complement.
    - Calculate the skew.
    - Find approximate kmers based on hamming distance.
'''

# Finding kmers

def kmer_count(seq, kmer):
    '''Counts occurences of a kmer in a sequence.'''
    count = 0
    for i in range(len(seq)-len(kmer)+1):
        if seq[i:i+len(kmer)] == kmer:
            count += 1
    return count

def frequency_map(seq, k):
    '''Finds the most frequent kmers of length k in a sequence.'''
    freq = {}
    n = len(seq)
    for i in range(n-k+1):  # first for loop selects each possible k-mer
        kmer = seq[i:i+k]
        freq[kmer] = 0
        for i in range(n-k+1):  # second loop searches Text for each k-mer
            if seq[i:i+k] == kmer:
                freq[kmer] = freq[kmer]+1
    return freq

def frequent_kmers(seq, k):
    '''Identifies the most frequent kmer in a sequence.'''
    kmers = []
    freq = frequency_map(seq, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            kmers.append(key)
    return kmers

def kmer_matching(seq, kmer):
    '''Identifies positions of a kmer in a sequence.'''
    positions = []
    n = len(seq)
    k = len(kmer)
    for i in range(n-k+1):
        if kmer == seq[i:i+k]:
            positions += [i]
    return positions


# Reverse complement
    
def reverse_complement(seq):
    '''Generates the reverse complement of a sequence'''
    seq = reverse(seq)
    seq = complement(seq)
    return seq

def reverse(seq):
    '''Reverses a sequence.'''
    rev = ''
    for base in seq:
        rev = base + rev
    return  rev

def complement(seq):
    '''Generates the complement of a sequence.'''
    comp = ''
    pairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    for base in seq:
        comp += pairs.get(base)
    return comp


# Skew
    
def symbol_array(genome, symbol):
    '''Calculates the relative proportions of a symbol in opposing halves
    of the genome.
    '''
    
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n//2]

    # look at the first half of genome to compute first array value
    array[0] = kmer_count(symbol, genome[0:n//2])

    for i in range(1, n):
        # set the current array value equal to the previous array value
        array[i] = array[i-1]

        # current array value can differ from previous value by at most 1
        if extended_genome[i-1] == symbol:
            array[i] = array[i]-1
        if extended_genome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def skew_array(genome):
    '''Calculates skew at every point in the genome.'''
    skew = [0]
    score = {'A':0, 'T':0, 'G':+1, 'C':-1}
    for i in range(0, len(genome)):
        skew.append(skew[i] + score[genome[i]])
    return skew

def minimum_skew(genome):
    '''Finds the position of minimum skew in the genome.'''
    positions = []
    array = skew_array(genome)
    m = min(array)
    for i in range(len(array)):
        if array[i] == m:
            positions.append(i)
    return positions


# kmer matching
    
def hamming_distance(p, q):
    '''Calculates the hamming distance between two sequences.'''
    ham = 0
    for i in range(len(p)):
        if q[i] != p[i]:
            ham += 1
    return ham

def approximate_kmer_matching(seq, kmer, ham):
    '''Finds positions of kmers within hamming distace.'''
    positions = []
    n = len(seq)
    k = len(kmer)
    for i in range(n-k+1):
        if hamming_distance(kmer, seq[i:i+k]) <= ham:
            positions.append(i)
    return positions

def approximate_kmer_count(kmer, seq, ham):
    '''Counts kmers in a sequence within hamming distance.'''
    count = 0
    for i in range(len(seq)-len(kmer)+1):
        if hamming_distance(kmer, seq[i:i+len(kmer)]) <=ham:
            count += 1
    return count