#Our MPRA has UMIs. UMIsperoligo takes a sam file of alignments, the fastq
#file containing the UMIs (on the reverse read) and gives a table of the following format:
#  oligo ID      number_of_reads     number_of_unique_UMIs

#So here we want to take that file (actually many of those files, one per condition per replicate)
#and calculate nucleotide-level enrichments. If an oligo covers a nucleotide, that nucleotide will
#be assigned one count. Then the idea is to compare depth-normalized counts for each nucleotide 
#across conditions.

#python3
import sys
import os
import numpy as np
from math import log2
import pandas as pd

def getntcoverage(umicounts, stepsize, oligosize, smoothcoverage):
	#Here we are going to get nt-level coverage. Oligo names are of the form
	#ENSMUSG00000025352.256|Gdf11
	#This tells you which UTR it belongs to (Gdf11) and where it is 
	#(it's number 256 in the UTR, counting 5' to 3')
	#If we know the step size, we can tell how far away this oligo is from the beginning
	#of the UTR. Further, if we know the oligo size, we know the nts it is covering.
	ntcoverage = {} #{utr : {nt : [normalized counts of oligos that cover that nt]}}  nt is 0-based
	medianntcoverage = {} #{utr : {nt : [median normalized counts of oligos that cover this nt]}}
	stepsize = int(stepsize)
	oligosize = int(oligosize)

	#Get total number of unique UMIs for depth normalization purposes
	with open(umicounts, 'r') as infh:
		counts = []
		for line in infh:
			uniqueumis = int(line.split('\t')[2])
			counts.append(uniqueumis)

		totalcounts = sum(counts) / 1e6


	with open(umicounts, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			utr = line[0].split('|')[1]
			if utr not in ntcoverage:
				ntcoverage[utr] = {}
			oligopos = int(line[0].split('.')[1].split('|')[0]) #this is 1-based
			oligostart = (oligopos - 1) * stepsize
			ntscovered = list(range(oligostart, oligostart + oligosize + 1))
			umicount = int(line[2])
			normalizedumicount = umicount / totalcounts
			for nt in ntscovered:
				if nt not in ntcoverage[utr]:
					ntcoverage[utr][nt] = [normalizedumicount]
				else:
					ntcoverage[utr][nt].append(normalizedumicount)

	for utr in ntcoverage:
		medianntcoverage[utr] = {}
		for nt in ntcoverage[utr]:
			mediancov = round(np.median(ntcoverage[utr][nt]), 3)
			medianntcoverage[utr][nt] = mediancov

	
	if smoothcoverage == False:
		return medianntcoverage
	elif smoothcoverage == True:
		#mean of window centered on nt
		smoothntcoverage = {} #{utr : {nt : smoothedcoverage}}
		windowsize = 11 #nt plus 5 nt on either side

		for utr in medianntcoverage:
			smoothntcoverage[utr] = {}
			maxpos = max(list(medianntcoverage[utr].keys()))
			for i in range(maxpos):
				if i < 5 or i > maxpos - 5: #window will be too small
					try:
						smoothntcoverage[utr][i] = medianntcoverage[utr][i]
					except KeyError:
						pass
				else:
					window = range(i - ((windowsize - 1) / 2), i + ((windowsize + 1) / 2))
					windowmedians = []
					for j in window:
						try:
							windowmedians.append(medianntcoverage[utr][j])
						except KeyError:
							pass
					smoothvalue = np.mean(windowmedians)
					smoothntcoverage[utr][i] = smoothvalue


		return smoothntcoverage

def collatereplicates(sampconds):
	#Given a file that tells how replicates are arranged, collate a 
	#pandas df of oligo umi counts
	samps = {}
	samps['conditionA'] = []
	samps['conditionB'] = []
	with open(sampconds, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'conditionA':
				continue
			samps['conditionA'].append(line[0])
			samps['conditionB'].append(line[1])

	condAcoverage = []
	condBcoverage = []
	for samp in samps['conditionA']:
		print('Calculating coverage for {0}...'.format(samp))
		coverage = getntcoverage(samp, 4, 260, True) #{utr : {nt : coverage}}
		condAcoverage.append(coverage)
	for samp in samps['conditionB']:
		print('Calculating coverage for {0}...'.format(samp))
		coverage = getntcoverage(samp, 4, 260, True)
		condBcoverage.append(coverage)

	#calculate log2FC at each nt
	log2fcdict = {} #{utr : {nt : log2FC (B/A)}}
	#define a reference coverage dict that we will use to loop through
	refdict = condAcoverage[0]
	for utr in refdict:
		log2fcdict[utr] = {}
		for nt in refdict[utr]:
			condAntcoverages = []
			condBntcoverages = []
			for samp in condAcoverage:
				try:
					coverage = samp[utr][nt]
					condAntcoverages.append(coverage)
				except KeyError:
					print(utr, nt)
					pass
			for samp in condBcoverage:
				try:
					coverage = samp[utr][nt]
					condBntcoverages.append(coverage)
				except KeyError:
					print(utr, nt)
					pass
			condAmeancoverage = np.mean(condAntcoverages)
			condBmeancoverage = np.mean(condBntcoverages)
			#add pseudocount
			pc = 1e-3
			condAmeancoverage = condAmeancoverage + pc
			condBmeancoverage = condBmeancoverage + pc
			log2FC = log2(condBmeancoverage / condAmeancoverage)
			log2fcdict[utr][nt] = log2FC

	#print(log2fcdict)
	
	#calculate pvalue for each nt


#Take 'sampconds' list of UMI count files
#collate ntcoverage dicts
#Take log2FC across condition of medians across replicates
#Repeat with 1000 shuffles
#Calculate empirical p value

#To define "regions":
#start at sig nt, allow max gap of <gap> nonsig nt
collatereplicates(sys.argv[1])


'''
x = getntcoverage(sys.argv[1], 4, 260, True)

with open('n2a20.txt', 'w') as outfh:
	outfh.write(('\t').join(['utr', 'nt', 'coverage']) + '\n')
	for utr in x:
		for nt in x[utr]:
			cov = str(x[utr][nt])
			outfh.write(('\t').join([utr, str(nt), cov]) + '\n')
'''