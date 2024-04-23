import random
from Bio.Seq import Seq
import numpy as np
from Bio import Entrez

Entrez.email = "ghassan.abboud@epfl.ch"
import random

def generate_weighted_sequence_variant(sequence, weights=[0.25, 0.25, 0.25, 0.25]):
    """
    Flips an initial sequence as it could happen when using the MinIon tool

    Arguments:
    sequence (Seq): sequence of interest
    weights(list of 4 float): probability of each sequence version (initial sequence, reverse complement, complement, reverse initial)

    Returns:
        selected_variant(Seq): randomly reversed sequence
        variant_description(str): string indicating which variant is selected
    """
    variants = ["initial sequence", "reverse complement", "complement", "reverse initial"]
    variant_sequences = [sequence, sequence.reverse_complement(), sequence.complement(), sequence.reverse_complement().complement()]
    
    selected_variant_index = random.choices(range(len(variants)), weights=weights)[0]
    selected_variant = variant_sequences[selected_variant_index]
    variant_description = variants[selected_variant_index]

    return selected_variant, variant_description


import random

def break_sequence_with_probability(sequence, break_prob_function):
    """
    Simulates breakages in the sequence, which could happen using the MinIon tool

    Arguments:
        sequence(Seq): sequence of interest
        break_prob_funtion(function): probability function used as weight for the breakages

    Returns:
        final_sequence (Seq): final sequence after breakages
        break_info (dict): dictionary containing information about the breakages
    """
    broken_sequence = []
    nbreaks = 2  # Max 2 breakages
    break_info = {'number_of_breaks': 0, 'part_taken': ''}
    
    for i, base in enumerate(sequence):
        # Computes the breakage probability
        break_prob = break_prob_function(i, len(sequence)) / 10

        # Checks if sequence breaks at this iteration
        if random.random() <= break_prob:
            if nbreaks > 0:
                broken_sequence.append('N')
                nbreaks -= 1   # Breakage marker is 'N'
                break_info['number_of_breaks'] += 1
        else:
            broken_sequence.append(base)

    if broken_sequence.count('N') == 2:  # If 2 breakages, take the middle part
        start_index = broken_sequence.index('N')
        end_index = len(broken_sequence) - broken_sequence[::-1].index('N') - 1
        final_sequence = broken_sequence[start_index + 1:end_index]
        break_info['part_taken'] = 'middle'
    elif broken_sequence.count('N') == 1:  # If 1 breakage, take the biggest part
        n_index = broken_sequence.index('N')
        if n_index < len(broken_sequence) - n_index - 1:
            final_sequence = broken_sequence[n_index + 1:]
            break_info['part_taken'] = 'end'
        else:
            final_sequence = broken_sequence[:n_index]
            break_info['part_taken'] = 'start'
    else:  # If 0 breakages, take full sequence
        final_sequence = broken_sequence
        break_info['part_taken'] = 'full sequence'

    return final_sequence, break_info


def break_prob_function(position, sequence_length):
    """
    Example of probability function usable for the previous function

    Arguments:
        position(int): position in sequence
        sequence_length(int)
    Returns:
        probability at a specific position in sequence(float)
    """
    max_prob = 0.5  # Max probability at start and end of sequence
    return max_prob * (position / sequence_length)

def mutation(base):
    """
    Simulates substitution mutation

    Arguments:
        base(char)

    Returns:
        muted base (char)
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
    
    Arguments:
        sequence(Seq): sequence of interest
        mutation_probability(float)
        mutation_mean(int): mean of normal law associated with mutation
        mutation_sd(int): standard deviation of normal law associated with mutation
        basic_mean(int): mean of normal law associated with no-change bases
        mutation_mean(int): standard deviation of normal law associated with no-change bases
    
    Returns:
        quality score of sequence (in base 64 characters -> ASCII)
        nb_mutations(int): number of mutations
        mutations_positions(array of ints): positions of mutations
    """
    quality_scores = []
    base64_table = "!°#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    nb_mutations=0
    mutations_positions= []
    for i, base in enumerate(sequence):
        # Computes score for each base depending on mutation
        if random.random() < mutation_probability:
            # If mutation needed
            base = mutation(base)
            nb_mutations+=1
            mutations_positions.append(i+1)
            quality_score = np.random.normal(mutation_mean, mutation_sd)
        else:
            # If no mutation needed
            quality_score = np.random.normal(basic_mean, basic_sd)
        # Limits the score between 0 and 93 (ASCII go from 33 to 126)
        quality_score = max(min(quality_score, 93), 0)
        base64_score = base64_table[int(round(quality_score))]
        quality_scores.append(base64_score)
    return ''.join(quality_scores), nb_mutations, mutations_positions

def download_sequence(species, gene_name, dst, start_length=None, stop_length= None, id = None, permissive_search = True):
    """
    Download sequence from GenBank through the Entrez database. 

    Parameters:
        species(str): name of species
        gene_name(str): name of gene
        dst(str,Path-like): destination file path
        start_length(int): minimum length of sequence
        stop_length(int): maximum length of sequence
        id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
        permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
    """
    
    if id == None:
        search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
        print(search_term)
        if start_length!= None or stop_length!= None:
                search_term += f" {start_length}:{stop_length}[Sequence Length]"
        handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
        search_result = Entrez.read(handle)
        handle.close()
        id = search_result["IdList"]
        n=0
        for i in id:
            n+=1
        if n==0 and permissive_search:
            search_term = f"{gene_name} {species} {start_length}:{stop_length}[Sequence Length]"
            handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
            search_result = Entrez.read(handle)
            handle.close()
            id = search_result["IdList"]
    n=0
    for i in id:
        n+=1
    if n==1:
        handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
        sequence = handle.read()
        handle.close()
        with open(dst, mode="a") as writer:
            writer.write(sequence)