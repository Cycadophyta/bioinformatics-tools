#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 10:04:47 2020

@author: pi
"""

'''
Motifs.py contains functions to identify motifs in parallel DNA sequences.

From a set of motifs, these functions can:
    - Generate a profile 
    - Generate a consensus sequence
    - Score the motifs
    - Calculte the probability of a kmer

From multiple sequences, these functions can search for the most probable 
motifs using three different methods:
    - Greedy motif search
    - Randomised motif search
    - Gibbs Sampling
'''

# Imports
import random


# Scoring motifs

def count(motifs):
    '''Counts bases in each position across motifs.'''
    k = len(motifs[0])
    count = {base:[0]*k for base in "ACGT"}
    for row in motifs:
        for idx, char in enumerate(row):
            count[char][idx] += 1
    return count

def profile(motifs):
    '''Generates a profile of motifs based on base counts.'''
    t = len(motifs)
    profile = count(motifs) 
    for key, val in profile.items():
        val[:] = [i / t for i in val]
    return profile

def consensus(motifs):
    '''Generates a consensus sequence based on base counts.'''
    k = len(motifs[0])
    counts = count(motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequent = ""
        for symbol in "ACGT":
            if counts[symbol][j] > m:
                m = counts[symbol][j]
                frequent = symbol
        consensus += frequent
    return consensus

def score(motifs):
    '''Scores motifs based on similarity to the consensus.'''
    con = consensus(motifs)
    score = 0
    for motif in motifs:
        for i in range(len(motif)):
            if motif[i] != con[i]:
                score += 1
    return score

def pr(kmer, profile):
    '''Calculates the probability of a kmer from a profile.'''
    pr = 1
    for i, base in enumerate(kmer):
        pr *= profile[base][i] 
    return pr

def probable_kmer(seq, k, profile):
    '''Returns the most probable kmer given a profile.'''
    score = -1  # score starts at -1 to allow for pr = 0
    for i in range(len(seq)-k+1):
        kmer = seq[i: i+k]
        kmer_pr = pr(kmer, profile)
        if kmer_pr > score:
            most = kmer
            score = kmer_pr
    return most


# Greedy motif search

def greedy_motif_search(dna, k, t):
    '''Finds the highest scoring combination of motifs.'''
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
    for i in range(len(dna[0])-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range (1, t):
            motifs.append(probable_kmer(dna[j], k, profile(motifs)))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


# Pseudocounts
    
def pseudocount(motifs):
    '''Counts bases in each position of motifs with pseudocounts.'''
    k = len(motifs[0])
    count = {base:[1]*k for base in "ACGT"}
    for row in motifs:
        for idx, base in enumerate(row):
            count[base][idx] += 1
    return count

def pseudoprofile(motifs):
    '''Calculates proportion of base counts in motifs with pseudocounts.'''
    profile = pseudocount(motifs)
    total = len(motifs) + 4
    for val in profile.values():
        val[:] = [i / total for i in val]
    return profile


# Greedy motif search with pseudocounts
    
def pseudo_greedy_motif_search(dna, k, t):
    '''Finds the highest scoring combination of motifs with pseudocounts.'''
    best_motifs = []
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])
    for i in range(len(dna[0])-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range (1, t):
            motifs.append(probable_kmer(dna[j], k, pseudoprofile(motifs)))
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


# Randomised motif search
    
def profile_motifs(profile, dna):
    '''Generates most probable motifs based on a profile.'''
    return[probable_kmer(seq, len(profile['A']), profile) for seq in dna]
    
def random_motifs(dna, k):
    '''Generates random motifs from given sequences.'''
    motifs = []
    for seq in dna:
        i = random.randint(0, len(seq)-k)
        motifs.append(seq[i:i+k])
    return motifs

def randomised_motif_search(dna, k):
    '''Finds the best set of motifs in dna.'''
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    while True:
        profile = pseudoprofile(motifs)
        motifs = profile_motifs(profile, dna)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs 
        
def best_randomised_motifs(dna, k, N):
    '''Repeats randomised_motif_search N times and 
    returns the best set of motifs.
    '''
    
    best_motifs = []
    best_score = float('inf')
    for n in range(N):
        motifs = randomised_motif_search(dna, k)
        if score(motifs) < best_score:
            best_score = score(motifs)
            best_motifs = motifs
    return best_motifs


# Gibbs sampling
    
def normalise(p):
    '''Normalise probabilities to add up to 1.'''
    return {key: p[key] / sum(p.values()) for key in p}
        
def weighted_die(p):
    '''Chooses a motif randomly, but respecting probabilities.'''
    return random.choices(
        list(p.keys()), 
        weights=list(p.values())
        )[0]
    
def profile_generated_kmer(seq, profile, k):
    '''Returns a random kmer from a sequence based on a profile.'''
    n = len(seq)
    p = {}
    for i in range(0, n-k+1):
        p[seq[i:i+k]] = pr(seq[i:i+k], profile)
    return weighted_die(normalise(p))

def gibbs_sampler(dna, k, t, N):
    '''Finds the best motifs based on Gibbs sampling.'''
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    for i in range(N):
        j = random.randint(0, t-1)
        profile = pseudoprofile(motifs[:j] + motifs[j+1:])
        motifs[j] = profile_generated_kmer(dna[j], profile, k)
        if score(motifs) > score(best_motifs):
            best_motifs = motifs
    return best_motifs

def best_gibbs_motifs(dna, k, t, N, repeats):
    '''Returns the best motifs from repeats of gibbs_sampler.'''
    best_motifs = []
    best_score = float('inf')
    for i in range(repeats):
        motifs = gibbs_sampler(dna, k, t, N)
        if score(motifs) < best_score:
            best_score = score(motifs)
            best_motifs = motifs
    return best_motifs