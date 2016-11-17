#!/usr/bin/env python

#script used to create ts and tn counts for SnIPRE
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import craig

dict = {'A': ('C','G','T' ), 'G': ('C','A','T'), 'C': ('T','A','G'), 'T': ('A','C','G')}

ref = {}
# create ref for the codon table for ts and tn
for codons in craig.ex.aminoAcidMap:
	ts = 0
	tn = 0
	base_aa = craig.ex.aminoAcidMap[codons]
	sites =  list(codons)
	for i in range(0,3):
		n = 0
		s = 0
		for j in dict[sites[i]]:
			if i == 0:
				new_codon = j + sites[1] + sites[2]	
				new_aa = craig.ex.aminoAcidMap[new_codon]	
				if new_aa != base_aa:
					n += 1
				else:
					s +=1
			if i == 1:
				new_codon = sites[0] + j  + sites[2]	
				new_aa = craig.ex.aminoAcidMap[new_codon]	
				if new_aa != base_aa:
					n += 1
				else:
					s +=1
			if i == 2:
				new_codon = sites[0] + sites[1]	+ j
				new_aa = craig.ex.aminoAcidMap[new_codon]	
				if new_aa != base_aa:
					n += 1
				else:
					s +=1
		if n > 0:
			tn += 1
		if s > 0:
			ts += 1
	
	#print "%s\t%s\t%s" % (codons, tn,ts)
	ref[codons]= (tn,ts)

file = sys.argv[1]
fasta_sequences = SeqIO.parse(open(file), 'fasta')
for fasta in fasta_sequences:
	ts = 0
	tn = 0
	name,description, sequence = fasta.id,fasta.description, str(fasta.seq)
	sequence = sequence.upper()
	codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
	if len(sequence) % 3 ==0:
		for i in codons:
			if i in craig.ex.aminoAcidMap:
				tmp_tn = 0
				tmp_ts = 0
				if i not in craig.ex.stopCodons:
					tmp_tn,tmp_ts = ref[i]	
					tn += tmp_tn
					ts += tmp_ts
	print "%s\t%s\t%s" % (name,tn,ts)
