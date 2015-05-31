#!/usr/bin/env python
codon2aa = {
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
'GGG': 'G', 'TAG': '*', 'TAA': '*', 'TGA':'*'}
aa2fullname = {
'F' : 'Phenylalanine',
'L' : 'Leucine',
'I' : 'Isoleucine',
'M' : 'Methionine',
'V' : 'Valine',
'S' : 'Serine',
'P' : 'Proline',
'T' : 'Threonine',
'A' : 'Alamine',
'Y' : 'Tyrosine',
'H' : 'Histidine',
'Q' : 'Glutamine',
'N' : 'Asparagine',
'K' : 'Lysine',
'D' : 'Aspartic acid',
'E' : 'Glutamic acid',
'C' : 'Cysteine',
'W' : 'Tryptophan',
'R' : 'Arginine',
'G' : 'Glycine',
'*': 'Stop Codon'}

def generate_codon_list():
    codon_list = []
    base_list = ['T', 'C', 'G', 'A']
    for b0 in base_list:
        for b1 in base_list:
            for b2 in base_list:
                codon = b0+b1+b2
                codon_list.append(codon)
    return codon_list

def generate_aa_list():
    codon_list = generate_codon_list()
    aa_list = []
    for codon in codon_list:
        aa = codon2aa[codon]
        if aa not in aa_list:
            aa_list.append(aa)
    return aa_list
