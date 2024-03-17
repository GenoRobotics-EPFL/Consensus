import random
from Bio.Seq import Seq
import numpy as np

def generate_weighted_sequence_variant(sequence, weights=[0.25, 0.25, 0.25, 0.25]):
    """
    flips an initial sequence as it could happen when using the MinIon tool
    :param sequence: sequence of interest [Seq]
    :param weights: probability of each sequence version (initial sequence, reverse complement, complement, reverse initial) [list of 4 [float]]
    :return: randomly reversed sequence [Seq]
    """
    #4 different outputs : 3'->5'/5'->3' and fwd/rev
    variants = [sequence, sequence.reverse_complement(), sequence.complement(), sequence.reverse_complement().complement()]
    # Random weighted choice
    selected_variant = random.choices(variants, weights=weights)[0]

    return selected_variant

def break_sequence_with_probability(sequence, break_prob_function):
    """
    simulates breakages in the sequence, which could happen using the MinIon tool
    :param sequence: sequence of interest [Seq]
    :param break_prob_funtion: probability function used as weight for the breakages [function]
    :return: final sequence after breakages
    """
    broken_sequence = []
    nbreaks=2 #Max 2 breakages 
    for i, base in enumerate(sequence):
        # Computes the breakage probability
        break_prob = break_prob_function(i, len(sequence)) / 10

        #Checks if sequence breaks at this iteration
        if random.random() <= break_prob:
            if nbreaks > 0 :
                broken_sequence.append('N')
                nbreaks-=1   #Breakage marker is 'N'
        else:
            broken_sequence.append(base)
    if broken_sequence.count('N') == 2:  # If 2 breakages, take the middle part
        start_index = broken_sequence.index('N')
        end_index = len(broken_sequence) - broken_sequence[::-1].index('N') - 1
        final_sequence = broken_sequence[start_index + 1:end_index]
    elif broken_sequence.count('N') == 1:  # If 1 breakage, take the biggest part
        n_index = broken_sequence.index('N')
        if n_index < len(broken_sequence) - n_index - 1:
            final_sequence = broken_sequence[n_index + 1:]
        else:
            final_sequence = broken_sequence[:n_index]
    else:  # If 0 breakages, take full sequence
        final_sequence = broken_sequence

    return Seq(''.join(final_sequence)) 

def break_prob_function(position, sequence_length):
    """
    Example of probability function usable for the previous function
    :param position: position in sequence [int]
    :param sequence_length: [int]
    :return: probability at a specific position in sequence [float]
    """
    max_prob = 0.5  # Max probability at start and end of sequence
    return max_prob * (position / sequence_length)

def mutation(base):
    """
    Simulates substitution mutation
    :param base: [char]
    :return: muted base [char]
    """
    # Bases list
    bases = ['A', 'T', 'C', 'G']

    # Randomly selected new base
    new_base = base
    while new_base == base:
        new_base = random.choice(bases)

    return new_base


def assign_quality_scores(sequence, mutation_probability=0.1,mutation_mean = 16,mutation_sd = 3,basic_mean = 48,basic_sd = 3):
    """
    Assigns fake scores for each base, following a different normal law if base is mutated or not
    :param sequence: sequence of interest [Seq]
    :param mutation_probability: [float]
    :param mutation_mean: mean of normal law associated with mutation [int]
    :param mutation_sd: standard deviation of normal law associated with mutation [int]
    :param mutation_mean: mean of normal law associated with no-change bases [int]
    :param mutation_mean: standard deviation of normal law associated with no-change bases [int]
    :return: quality score of sequence (in base 64 characters -> ASCII)
    """
    quality_scores = []
    base64_table = "!°#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    for base in sequence:
        # Computes score for each base depending on mutation
        if random.random() < mutation_probability:
            # If mutation needed
            base = mutation(base)
            quality_score = np.random.normal(mutation_mean, mutation_sd)
        else:
            # If no mutation needed
            quality_score = np.random.normal(basic_mean, basic_sd)
        # Limits the score between 0 and 93 (ASCII go from 33 to 126)
        quality_score = max(min(quality_score, 93), 0)
        base64_score = base64_table[int(round(quality_score))]
        quality_scores.append(base64_score)
    return ''.join(quality_scores)