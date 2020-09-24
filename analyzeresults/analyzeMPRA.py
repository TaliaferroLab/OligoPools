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
import numpy as np

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
			#normalizedumicount = umicount / totalcounts
			for nt in ntscovered:
				if nt not in ntcoverage[utr]:
					ntcoverage[utr][nt] = [umicount]
				else:
					ntcoverage[utr][nt].append(umicount)

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
					window = range(i - int(((windowsize - 1) / 2)), i + int(((windowsize + 1) / 2)))
					windowmedians = []
					for j in window:
						try:
							windowmedians.append(medianntcoverage[utr][j])
						except KeyError:
							pass
					smoothvalue = np.mean(windowmedians)
					smoothntcoverage[utr][i] = round(smoothvalue, 3)


		return smoothntcoverage

def collatereplicates(sampconds):
	#Given a file that tells how replicates are arranged, collate a 
	#pandas df of oligo umi counts
	samps = {} #{condition : [sample1, sample2, ...]}
	samps['conditionA'] = []
	samps['conditionB'] = []
	with open(sampconds, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'conditionA':
				continue
			samps['conditionA'].append(line[0])
			samps['conditionB'].append(line[1])

	allsamps = samps['conditionA'] + samps['conditionB']
	sampnames = [os.path.basename(samp) for samp in allsamps]
	coverages = [] #list of coverage dictionaries, samples in same order as allsamps

	#For each sample, make a one column dataframe where the rows are UTR-nt positions
	allsampdfs = []
	for idx, samp in enumerate(allsamps):
		print('Getting UMI counts from {0}...'.format(samp))
		umidict = {} #{utr : {oligonumber : umicount}}
		with open(samp, 'r') as infh:
			for line in infh:
				line = line.strip().split('\t')
				utr = line[0].split('|')[1]
				oligonumber = int(line[0].split('.')[1].split('|')[0])
				umicount = int(line[2])
				if utr not in umidict:
					umidict[utr] = {}
				umidict[utr][oligonumber] = umicount

		dfs = []
		for utr in umidict:
			df = pd.DataFrame.from_dict([umidict[utr]], orient = 'columns')
			cnames = list(df.columns.values)
			cnames = [utr + '_' + str(c) for c in cnames]
			df.columns = cnames
			df = df.transpose()
			df.columns = [sampnames[idx]]
			dfs.append(df)
		df = pd.concat(dfs, axis = 0)
		allsampdfs.append(df)

	#Perform outer join to bring all sample-level dfs together
	bigdf = pd.concat(allsampdfs, axis = 1)
	#There will be some NAs following this outer join, so replace the NAs with 0
	bigdf = bigdf.fillna(value = 0.0)

	return bigdf


def normalizecounts(countdf):
	print('Quantile normalizing...')
	#Given a dataframe of umi counts per sample, perform quantile normalization
	#See https://ai.recodeminds.com/news/computing-quantile-normalization-in-python/

	#Sort each column
	countdf_sorted = pd.DataFrame(np.sort(countdf.values, axis = 0), index = countdf.index, columns = countdf.columns)
	#Compute row means
	countdf_mean = countdf_sorted.mean(axis = 1)
	countdf_mean.index = np.arange(1, len(countdf_mean) + 1)
	#Use average values for each sample in the original order
	countdf_qn = countdf.rank(method = 'min').stack().astype(int).map(countdf_mean).unstack()

	return countdf_qn

def getntcoverage(countdf, stepsize, oligosize):
	#Here we are going to get nt-level coverage. Oligo names are of the form
	#ENSMUSG00000025352.256|Gdf11
	#This tells you which UTR it belongs to (Gdf11) and where it is 
	#(it's number 256 in the UTR, counting 5' to 3')
	#If we know the step size, we can tell how far away this oligo is from the beginning
	#of the UTR. Further, if we know the oligo size, we know the nts it is covering.

	#At the end of this, we want a df that has the median number of (normalized) oligo
	#counts across all oligos that cover a nt

	countdfindicies = list(countdf.index.values) #utr_oligonumber
	countdfcolnames = list(countdf.columns.values) #sample name
	allsampdfs = []

	for colname in countdfcolnames:
		print('Getting nt level counts for {0}...'.format(colname))
		ntcoverage = {} #{utr : {nt : [(normalized) counts of oligos that cover that nt]}}  nt is 0-based
		medianntcoverage = {} #{utr : {nt : median counts of oligos that cover this nt}}
		counts = countdf[colname].tolist()
		#indicies (oligonames) are utr_oligonumber (starts at 1)
		for idx, oligo in enumerate(countdfindicies):
			oligocounts = counts[idx]
			utr = oligo.split('_')[0]
			if utr not in ntcoverage:
				ntcoverage[utr] = {}
			oligonumber = int(oligo.split('_')[1])
			oligostart = (oligonumber - 1) * stepsize
			ntscovered = list(range(oligostart, oligostart + oligosize + 1))
			for nt in ntscovered:
				if nt not in ntcoverage[utr]:
					ntcoverage[utr][nt] = [oligocounts]
				else:
					ntcoverage[utr][nt].append(oligocounts)

		for utr in ntcoverage:
			medianntcoverage[utr] = {}
			for nt in ntcoverage[utr]:
				mediancov = round(np.median(ntcoverage[utr][nt]), 3)
				medianntcoverage[utr][nt] = mediancov

		#Turn this back into a df, one per utr
		utrdfs = []
		for utr in medianntcoverage:
			df = pd.DataFrame.from_dict([medianntcoverage[utr]], orient = 'columns')
			cnames = list(df.columns.values)
			cnames = [utr + '_' + str(c) for c in cnames]
			df.columns = cnames
			df = df.transpose()
			df.columns = [colname]
			utrdfs.append(df)

		#Bring all utr dfs together for this sample
		df = pd.concat(utrdfs, axis = 0)
		allsampdfs.append(df)

	#Perform outer join to bring all sample-level dfs together
	bigdf = pd.concat(allsampdfs, axis = 1)
	#There will be some NAs following this outer join, so replace the NAs with 0
	bigdf = bigdf.fillna(value = 0.0)
	
	return bigdf


def smoothcounts(ntdf, windowsize):
	#Given a df of counts, smooth counts by sliding a window centered on nt and taking the mean
	#of counts in that window

	#first nt of a utr is 0

	#Window size must be odd so that it can be centered on one nt
	windowsize = int(windowsize)
	if windowsize % 2 == 0:
		print('ERROR: window size must be odd!')
		sys.exit()

	ntdfindicies = list(ntdf.index.values) #utr_ntnumber
	ntdfcolnames = list(ntdf.columns.values) #sample name
	allsampdfs = []

	for colname in ntdfcolnames: #for every sample
		print('Smoothing counts for {0}...'.format(colname))
		ntcounts = {} #{utr : {nt : counts}}
		smoothedntcounts = {} #{utr : {nt : smoothedcounts}}
		counts = ntdf[colname].tolist()
		for idx, ntid in enumerate(ntdfindicies):
			ntc = counts[idx]
			utr = ntid.split('_')[0]
			if utr not in ntcounts:
				ntcounts[utr] = {}
			ntnumber = int(ntid.split('_')[1])
			ntcounts[utr][ntnumber] = ntc

		#now smooth this dictionary of counts
		for utr in ntcounts:
			smoothedntcounts[utr] = {}
			maxpos = max(list(ntcounts[utr].keys()))
			for i in range(maxpos + 1): #for every nt in the utr
				windowstart = int(max([0, i - ((windowsize - 1) / 2)]))
				windowend = int(min([maxpos, i + ((windowsize - 1) / 2)]))

				window = range(windowstart, windowend + 1)
				windowcounts = []
				for j in window:
					windowcounts.append(ntcounts[utr][j])
				smoothvalue = np.mean(windowcounts)
				smoothedntcounts[utr][i] = round(smoothvalue, 3)

		#turn this back into a df, one per utr
		utrdfs = []
		for utr in smoothedntcounts:
			df = pd.DataFrame.from_dict([smoothedntcounts[utr]], orient = 'columns')
			cnames = list(df.columns.values)
			cnames = [utr + '_' + str(c) for c in cnames]
			df.columns = cnames
			df = df.transpose()
			df.columns = [colname]
			utrdfs.append(df)

		#Bring all utr dfs together for this sample
		df = pd.concat(utrdfs, axis = 0)
		allsampdfs.append(df)

	#Perform outer join to bring all sample-level dfs together
	bigdf = pd.concat(allsampdfs, axis = 1)

	return bigdf




#Take 'sampconds' list of UMI count files <done>
#collate ntcoverage dicts <done>
#Take log2FC across condition of medians across replicates
#Repeat with 1000 shuffles
#Calculate empirical p value

#To define "regions":
#start at sig nt, allow max gap of <gap> nonsig nt
umidf = collatereplicates(sys.argv[1])
umidf_qn = normalizecounts(umidf)
umidf.to_csv(path_or_buf = 'bigdf.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')
umidf_qn.to_csv(path_or_buf = 'bigdf_qn.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')
ntdf = getntcoverage(umidf_qn, 4, 260)
ntdf.to_csv(path_or_buf = 'bigdf_qn_nt.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')
smoothntdf = smoothcounts(ntdf, 11)
smoothntdf.to_csv(path_or_buf = 'bigdf_qn_nt_smooth.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')

'''
reorder steps
1. get umi counts
2. quantile normalize
3. get nuc counts
4. optionally, smooth
'''