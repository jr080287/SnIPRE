#!/usr/bin/env python

"""Import moduels"""
import sys
from Bio import AlignIO
from Bio.Seq import Seq
input_file = sys.argv[1]
outgroup_n = int(sys.argv[2])
ingroup_codon = {}
ingroup_aa = {}
ind_count = 0
alignment = AlignIO.read(open(input_file), 'fasta')
total_len = alignment.get_alignment_length()
individuals = len(alignment)

outgroup_aln = alignment[individuals-outgroup_n:individuals]
ingroup_aln = alignment[0:(individuals - outgroup_n)]

polysites = {}
for site in range(0, total_len):
    nucs = {}
    for record in alignment:
        if 'N' not in record.seq[site] and '-' not in record.seq[site]:
            nucs[record.seq[site]] = 1

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
        tmp_codons = {}
        codon_start = (codonid - 1) * 3
        codon_end = (codon_start + 3)
        codonid = sites
        for record in ingroup_aln:
            codon = str(record.seq[codon_start:codon_end])
            if codon not in tmp_codons and 'N' not in codon and '-' not in codon:
                tmp_codons[codon] = 1
        tmp_aa = {}
        for i in tmp_codons:
            aa = str(Seq(str(i)).translate())
            if aa not in tmp_aa:
                tmp_aa[aa] = 1
        if len(tmp_codons) > 1 and len(tmp_codons) < 3: #this is a polymorphic site within the ingroup
            if len(tmp_aa) > 1:
                pn += 1
            else:
                ps += 1
        else: #this is a fixed within the ingroup
            if outgroup_n < 2:
                codon = outgroup_aln.seq[codon_start:codon_end]
                if 'N' not in codon and '.' not in codon and '-' not in codon:
                    aa = str(Seq(str(codon)).translate())
                    if aa in tmp_aa:
                        ds += 1
                    else:
                        dn += 1
            else:
                tmp_codons_out = {}
                for record in outgroup_aln:
                    codon = str(record.seq[codon_start:codon_end])
                    if codon not in tmp_codons and 'N' not in codon and '-' not in codon:
                        tmp_codons_out[codon] = 1
                tmp_aa_out = {}
                for i in tmp_codons_out:
                    aa = str(Seq(str(i)).translate())
                    if aa not in tmp_aa_out:
                        tmp_aa_out[aa] = 1
                if len(tmp_codons_out) > 1: #this is polymorphic in the outgroup
                    if len(tmp_aa_out) > 1:
                        #make sure it is not triallic at the codon level bc of the outgroup
                        tmp = list(set(tmp_codons.keys()+tmp_codons_out.keys()))
                        if tmp > 1 and tmp < 3:
                            if len(list(set(tmp_aa.keys()+tmp_aa_out.keys()))) > 1:
                                pn += 1
                            else:
                                ps += 1
                else: #this is fixed within the outgroup
                    if (len(list(set(tmp_aa.keys() + tmp_aa_out.keys())))) > 1:
                        dn += 1
                    else:
                        ds += 1

print "%s\t%s\t%s\t%s" % (dn, ds, pn, ps)
