from prody import *
from Bio import pairwise2
from Bio import PDB
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from subprocess import call
import numpy as np
from joblib import Parallel, delayed
import multiprocessing as mp
import gc
import time

thold = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]

def coevol(line):
	result={}
	entry = line.split("\t")
	bound_pro = entry[0].split("_")[0]
	bound_pro_c1 = entry[0].split("_")[1].split(":")[0]
	bound_pro_c2 = entry[0].split("_")[1].split(":")[1]

	unbound_pro1 = entry[1].split("_")[0]
	unbound_pro1_c = entry[1].split("_")[1]

	unbound_pro2 = entry[2].split("_")[0]
	unbound_pro2_c = entry[2].split("_")[1][:-1]

	header = parsePDBHeader(bound_pro, 'polymers')
	bp_chid = []
	for polymer in header:
		bp_chid.append(polymer.chid)

	bp1_chid_idx = []
	for c in bound_pro_c1:
		if c in bp_chid:
			bp1_chid_idx.append(bp_chid.index(c))

	bp2_chid_idx = []
	for c in bound_pro_c2:
		if c in bp_chid:
			bp2_chid_idx.append(bp_chid.index(c))

	# header = parsePDBHeader(unbound_pro1, 'polymers')
	# up1_chid = []
	# for polymer in header:
	# 	up1_chid.append(polymer.chid)

	# up1_chid_idx = []
	# for c in unbound_pro1_c:
	# 	if c in up1_chid:
	# 		up1_chid_idx.append(up1_chid.index(c))

	# header = parsePDBHeader(unbound_pro2, 'polymers')
	# up2_chid = []
	# for polymer in header:
	# 	up2_chid.append(polymer.chid)

	# up2_chid_idx = []
	# for c in unbound_pro2_c:
	# 	if c in up2_chid:
	# 		up2_chid_idx.append(up2_chid.index(c))
	seqs = []
	pfam1 = []
	for i in range(len(bp1_chid_idx)):
		if i == 0:
			unip_raw=str(header[bp1_chid_idx[i]].dbrefs).split(" ")[1]
			try:
				unip = searchUniprotID(unip_raw)
				pfamid=searchPfam(unip).keys()[0]
				seq = pdb.getHierView()[bp_chid[bp1_chid_idx[i]]].getSequence()
				good = 1
			except IndexError:
				pdb = parsePDB(bound_pro)
				seq = pdb.getHierView()[bp_chid[bp1_chid_idx[i]]].getSequence()
				pfamid=searchPfam(seq).keys()[0]
				good = 0

			fetchPfamMSA(pfamid)
			pfam1.append(pfamid)
			raw_msa = parseMSA(pfamid + '_full.sth')
			if good == 1:
				refined_msa = refineMSA(raw_msa, label=unip)
			else:
				refined_msa = raw_msa
			total_msa = refined_msa
			total_seq = seq
			seqs.append(seq)
		else:
			unip_raw=str(header[bp1_chid_idx[i]].dbrefs).split(" ")[1]
			try:
				unip = searchUniprotID(unip_raw)
				pfamid=searchPfam(unip).keys()[0]
				seq = pdb.getHierView()[bp_chid[bp1_chid_idx[i]]].getSequence()
				good = 1
			except IndexError:
				pdb = parsePDB(bound_pro)
				seq = pdb.getHierView()[bp_chid[bp1_chid_idx[i]]].getSequence()
				pfamid=searchPfam(seq).keys()[0]
				good = 0
			fetchPfamMSA(pfamid)
			pfam1.append(pfamid)
			raw_msa = parseMSA(pfamid + '_full.sth')
			if good == 1:
				refined_msa = refineMSA(raw_msa, label=unip)
			else:
				refined_msa = raw_msa
			total_msa = mergeMSA(total_msa, refined_msa)
			total_seq = total_seq + seq
			seqs.append(seq)
	total_msa1 = total_msa
	total_seq1 = total_seq

	pfam2 = []
	for i in range(len(bp2_chid_idx)):
		if i == 0:
			unip_raw=str(header[bp2_chid_idx[i]].dbrefs).split(" ")[1]
			try:
				unip = searchUniprotID(unip_raw)
				pfamid=searchPfam(unip).keys()[0]
				seq = pdb.getHierView()[bp_chid[bp2_chid_idx[i]]].getSequence()
				good = 1
			except IndexError:
				pdb = parsePDB(bound_pro)
				seq = pdb.getHierView()[bp_chid[bp2_chid_idx[i]]].getSequence()
				pfamid=searchPfam(seq).keys()[0]
				good = 0

			fetchPfamMSA(pfamid)
			pfam2.append(pfamid)
			raw_msa = parseMSA(pfamid + '_full.sth')
			if good == 1:
				refined_msa = refineMSA(raw_msa, label=unip)
			else:
				refined_msa = raw_msa
			total_msa = refined_msa
			total_seq = seq
			seqs.append(seq)
		else:
			unip_raw=str(header[bp2_chid_idx[i]].dbrefs).split(" ")[1]
			try:
				unip = searchUniprotID(unip_raw)
				pfamid=searchPfam(unip).keys()[0]
				seq = pdb.getHierView()[bp_chid[bp2_chid_idx[i]]].getSequence()
				good = 1
			except IndexError:
				pdb = parsePDB(bound_pro)
				seq = pdb.getHierView()[bp_chid[bp2_chid_idx[i]]].getSequence()
				pfamid=searchPfam(seq).keys()[0]
				good = 0
			fetchPfamMSA(pfamid)
			pfam1.append(pfamid)
			raw_msa = parseMSA(pfamid + '_full.sth')
			if good == 1:
				refined_msa = refineMSA(raw_msa, label=unip)
			else:
				refined_msa = raw_msa
			total_msa = mergeMSA(total_msa, refined_msa)
			total_seq = total_seq + seq
			seqs.append(seq)
	total_msa2 = total_msa
	total_seq2 = total_seq

	mergedMSA = specMergeMSA(total_msa1, total_msa2)
	finalMergedMSA = refineMSA(mergedMSA, colocc=0.9)
	writeMSA(bound_pro + '.fasta', finalMergedMSA)

	merged_seq = total_seq1 + total_seq2
	g = open(bound_pro + '_one.fasta','w')
	g.write(">" + bound_pro + "\n")
	g.write(merged_seq)
	g.close()
	
	call(["clustalw -profile1=" + bound_pro + '.fasta -sequences -profile2=' + bound_pro + '_one.fasta'], shell=True)
	AlignIO.convert(bound_pro + "_one.aln","clustal",bound_pro + "_final.fasta","fasta")
	finalMSA = parseMSA(bound_pro + "_final.fasta")
	finalRefinedMSA = refineMSA(finalMSA, colocc=0.95)

	idx_real = finalRefinedMSA.getIndex(bound_pro)
	seq_from_alignment = str(finalRefinedMSA[idx_real])

	res_idx_str = []
	res_idx_ali = []
	for seq in seqs:
		alignment = pairwise2.align.globalms(seq, seq_from_alignment,5,-1,-10,-1)
		p1res = []
		p2res = []
		count = 0
		for i in range(0,len(alignment[0][0])):
				if alignment[0][0][i] == '-':
					count = count+1
				if alignment[0][1][i] == '-':
					count1 = count1 + 1
				if alignment[0][0][i] != '-' and alignment[0][1][i] != '-':
					p1res.append(i-count)
					p2res.append(i)
		res_idx_str.append(p1res)
		res_idx_ali.append(p2res)



	
	#call(["clustalw -infile=" + bound_pro + "_one.aln -outfile=" + bound_pro + "_final.fasta"], shell=True)

