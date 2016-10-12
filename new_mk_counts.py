#!/usr/bin/env python

import sys
from Bio.Seq import translate
from Bio import AlignIO
from Bio.Seq import Seq
input_file = sys.argv[1]
ingroup_codon = {}
ingroup_aa  ={}
ind_count = 0
alignment=AlignIO.read(open(input_file),'fasta')
total_len = alignment.get_alignment_length()
individuals = len(alignment)
outgroup_aln = alignment[-1]
ingroup_aln =  alignment[0:(individuals -1)]

polysites={}
for site in range(0,total_len):
    nucs = {}
    for record in alignment:
        if 'N' not in record.seq[site] and '-' not in record.seq[site]:
            nucs[record.seq[site]]=1

    if len(nucs) > 1:
        codonid = int((site/3))+1
        if codonid in polysites:
            polysites[codonid] += 1
        else:
            polysites[codonid] = 1

dn = 0
ds = 0
pn = 0
ps = 0
for sites in polysites:
    if polysites[sites] < 2:
        tmp_codons={}
        codon_start = (codonid - 1) * 3
        codon_end = (codon_start + 3 )
        codonid = sites
        for record in ingroup_aln:
            codon = str(record.seq[codon_start:codon_end])
            if codon not in tmp_codons and 'N' not in codon and '-' not in codon:
                tmp_codons[codon]=1
        tmp_aa = {}
        if len(tmp_codons) > 1: #this is a polymorphic site within the ingroup
            for i in tmp_codons:
                aa = str(Seq(str(i)).translate())
                if aa not in tmp_aa:
                    tmp_aa[aa]=1
            if len(tmp_aa) > 1:
                pn += 1
            else:
                ps += 1
        else: #this is a fixed within the ingroup
            codon = outgroup_aln.seq[codon_start:codon_end]
            if 'N' not in codon and '.' not in codon and '-' not in codon:
                aa = str(Seq(str(codon)).translate())
                if aa in tmp_aa:
                    ds += 1
                else:
                    dn += 1




print "%s\t%s\t%s\t%s" % (dn, ds, pn, ps)

