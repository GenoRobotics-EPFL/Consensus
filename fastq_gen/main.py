#main goal : generate "fake" fastq files from inital sequence (fasta file) with modulable parameters
#parameters that are subject to changes : sequence orientation (5'->3'/3'->5' and fwd/rev), sequence breaking, alignment scores, mutations
import os
import os.path as ospath
import utils as utl
from Bio.Seq import Seq
#Testing all coded functions with an example sequence and 10 iterations, the fastq file is generated in the outputs folder

sequence = Seq("ATCGGCTAGGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")
niter = 10
id = "@Sequence1"
if not ospath.exists("outputs/"):
    os.makedirs("outputs/")
with open(ospath.join(os.getcwd(), "outputs","output.fastq"), "a") as handle: # Operations on newly created or updated file
    for i in range(niter) :
        handle.write(f"@read{i}\n") #writes the sequence id
        seq = utl.generate_weighted_sequence_variant(sequence,[0.05,0.25,0.5,0.25])
        broken_seq = utl.break_sequence_with_probability(seq,utl.break_prob_function)
        quality_score = utl.assign_quality_scores(broken_seq,0.5,10,2,50,5)
        handle.write(str(broken_seq)+ '\n') #writes the sequence generated by transforming it into string
        handle.write("+\n") #writes the + separator
        handle.write(quality_score + '\n') #writes the quality score