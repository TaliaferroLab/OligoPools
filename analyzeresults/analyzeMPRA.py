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
import gffutils
import numpy as np
from math import log2
import pandas as pd
import numpy as np
from collections import OrderedDict
from statsmodels.robust.scale import mad as MAD
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from random import shuffle
from scipy.stats import percentileofscore
from scipy.stats.distributions import chi2
import warnings
import time

#go from nt number to genome coord, need {oligo : [list of exonic coordinates covered by oligos]}
#from that, make dictionary of {oligo : {utr nt(as defined by getntcoverage) : genome coord}}

def oligos2genome(gff):
	#We need some way to relate oligo positions, namely through their ID's (e.g. ENSMUSG00000025352.256|Gdf11)
	#to positions in the genome that they cover. One way to do this would be to, for each oligo,
	#make a list of coordinates that are tiled by the oligo: {utr : [list of exonic coordinates covered by utr's oligos]}
	#then, from that, make dictionary of {utr : {utr nt (as defined by getntcoverage (this goes from 1 to len(utr))) : genome coord}}

	#Gffs are annotations of oligo positions and were made by OligoPools_shortstep_260nt.py

	utrcoords = {} #{utr : [list of exonic coordinates covered by utr's oligos]}

	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	utrs = db.features_of_type('UTR')

	for utr in utrs:
		utrid = str(utr.attributes['gene_name'][0])
		utrexoniccoords = []
		for oligo in db.children(utr, featuretype = 'oligo', level = 1):
			#Does this oligo cross a splice junction (i.e. does it have junctionpiece children?)
			if len(list(db.children(oligo, featuretype = 'junctionpiece'))) > 0:
				#yes there is a junctionpiece
				for jp in db.children(oligo, featuretype = 'junctionpiece'):
					jpnt = list(range(jp.start, jp.end + 1))
					utrexoniccoords += jpnt
			elif len(list(db.children(oligo, featuretype = 'junctionpiece'))) == 0:
				#no there is not a junctionpiece
				oligont = list(range(oligo.start, oligo.end + 1))
				utrexoniccoords += oligont

		#there's gonna be a lot of redundancy in utrexoniccoords because the oligos overlap
		utrexoniccoords = list(set(utrexoniccoords))
		if utr.strand == '+':
			utrexoniccoords = list(sorted(utrexoniccoords))
		elif utr.strand == '-':
			utrexoniccoords = list(reversed(sorted(utrexoniccoords)))

		utrcoords[utrid] = utrexoniccoords

	#OK now we need to go from utrcoords to utr_oligocoords {utr : {utr nt (as defined by getntcoverage (this goes from 1 to len(utr))) : genome coord}}

	utr_oligocoords = {} #{utr : {utr nt (as defined by getntcoverage (this goes from 1 to len(utr))) : genome coord}}
	for utr in utrcoords:
		if utr != 'Gdf11':
			continue
		genomecoords = utrcoords[utr]
		oligocoords = list(range(len(genomecoords) + 1))
		oligocoords2genomecoords = dict(zip(oligocoords, genomecoords))
		utr_oligocoords[utr] = oligocoords2genomecoords

	
	#turn utr_oligocoords into a df and return the df
	df = pd.concat({k: pd.DataFrame.from_dict(v, 'index') for k, v in utr_oligocoords.items()}, axis = 0)
	#get indexes out as column names
	df.reset_index(level = 1, inplace = True)
	df.reset_index(level = 0, inplace = True)
	df.columns = ['Gene', 'utrcoord', 'genomecoord']

	return df


def collatereplicates(sampconds, conditionA, conditionB):
	#Given a file that tells how replicates are arranged, collate a 
	#pandas df of oligo umi counts
	samps = {} #{condition : [sample1, sample2, ...]}
	samps['conditionA'] = []
	samps['conditionB'] = []
	with open(sampconds, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'sample':
				continue
			if line[1] == conditionA:
				samps['conditionA'].append(line[0])
			elif line[1] == conditionB:
				samps['conditionB'].append(line[0])

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

	#Don't use smoothing here, 

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

def getlog2fc(countdf, sampconds, conditionA, conditionB):
	#Given a df of nucleotide-level counts and a file telling which samples 
	#belong to which condition, calculate mean counts for each condition
	#across replicates and then log2FC between conditions
	#log2FC are recorded as B/A

	#countdf is a df where indicies (rownames) are utr_nt and colnames are samples
	print('Calculating log2FC...')

	log2fcdict = {} #{utr : {nt : log2FC}}

	samps = {} #{condition : [sample1, sample2, ...]}
	samps['conditionA'] = []
	samps['conditionB'] = []
	with open(sampconds, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'sample':
				continue
			if line[1] == conditionA:
				samps['conditionA'].append(os.path.basename(line[0]))
			elif line[1] == conditionB:
				samps['conditionB'].append(os.path.basename(line[0]))

	for index, row in countdf.iterrows():
		utr = index.split('_')[0]
		nt = int(index.split('_')[1])
		if utr not in log2fcdict:
			log2fcdict[utr] = {}

		conditionAcounts = []
		conditionBcounts = []
		for sample in samps['conditionA']:
			counts = row[sample]
			conditionAcounts.append(counts)
		for sample in samps['conditionB']:
			counts = row[sample]
			conditionBcounts.append(counts)

		conditionAmeancounts = np.mean(conditionAcounts)
		conditionBmeancounts = np.mean(conditionBcounts)

		pc = 1 #pseudocount
		log2fc = log2((conditionBmeancounts + pc) / (conditionAmeancounts + pc))

		log2fcdict[utr][nt] = round(log2fc, 3)

	#turn this into a df, one per utr
	utrdfs = []
	for utr in log2fcdict:
		df = pd.DataFrame.from_dict([log2fcdict[utr]], orient = 'columns')
		cnames = list(df.columns.values)
		cnames = [utr + '_' + str(c) for c in cnames]
		df.columns = cnames
		df = df.transpose()
		df.columns = ['log2FC']
		utrdfs.append(df)

	#Bring together all utrdfs
	allutrdfs = pd.concat(utrdfs, axis = 0)
	allutrdfs.columns = ['log2FC']

	return allutrdfs


def getpvalues_MADstat(countdf, sampconds, conditionA, conditionB):
	#For each nt, calculate a stat (difference between medians / median absolute deviation)
	#median absolute deviation is the absolute median of deviations from the median

	#Uses log2-transformed values

	#shuffle sample labels, calculate stat again
	#repeat shuffle 1000 times
	#calculate empirical p value

	#alternatively, maybe use a LME model?

	#in countdf, indicies (rownames) are utr_nt and columns are samples

	samps = {} #{condition : [sample1, sample2, ...]}
	samps['conditionA'] = []
	samps['conditionB'] = []
	with open(sampconds, 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			if line[0] == 'sample':
				continue
			if line[1] == conditionA:
				samps['conditionA'].append(os.path.basename(line[0]))
			elif line[1] == conditionB:
				samps['conditionB'].append(os.path.basename(line[0]))

	pvaluedict = OrderedDict() #{utr_nt : p}
	condameandict = {} #{utr_nt : condAmean}
	condbmeandict = {} #{utr_nt : condAmean}
	condaMADdict = {} #{utr_nt : condAmad}
	condbMADdict = {} #{utr_nt : condBmad}
	currentutr = ''
	for index, row in countdf.iterrows():
		utr = index.split('_')[0]
		if utr != currentutr:
			print('Calculating pvalues for {0}...'.format(utr))
			currentutr = utr
		if utr != 'Gdf11':
			continue

		conditionAvalues = countdf.loc[index, samps['conditionA']].tolist()
		conditionBvalues = countdf.loc[index, samps['conditionB']].tolist()

		#log2transform
		pc = 0.01
		conditionAvalues = [log2(v + pc) for v in conditionAvalues]
		conditionBvalues = [log2(v + pc) for v in conditionBvalues]

		condAmean = np.mean(conditionAvalues)
		condBmean = np.mean(conditionBvalues)

		#by default, statsmodels.robust.scale.mad scales the MAD by dividing it by 0.674
		#which gives you approximately the standard deviation. Override by supplying c = 1

		condAmad = MAD(conditionAvalues, c = 1)
		condBmad = MAD(conditionBvalues, c = 1)

		a = (condAmad**2) / len(conditionAvalues)
		b = (condBmad**2) / len(conditionBvalues)
		bottom = (a + b)**0.5
		top = condBmean - condAmean
		stat = top / bottom

		allvalues = conditionAvalues + conditionBvalues
		nullstats = []
		for i in range(100):
			#shuffle values
			shuffle(allvalues)
			shufflecondAvalues = allvalues[:len(conditionAvalues)] #first n are new condA
			shufflecondBvalues = allvalues[len(conditionAvalues):] #everything else is condB

			shufflecondAmean = np.mean(shufflecondAvalues)
			shufflecondBmean = np.mean(shufflecondBvalues)
			
			shufflecondAmad = MAD(shufflecondAvalues, c = 1)
			shufflecondBmad = MAD(shufflecondBvalues, c = 1)

			a = (shufflecondAmad**2) / len(shufflecondAvalues)
			b = (shufflecondBmad**2) / len(shufflecondBvalues)
			bottom = (a + b)**0.5
			top = shufflecondBmean - shufflecondAmean
			nullstat = top / bottom
			nullstats.append(nullstat)

		
		percentile = percentileofscore(nullstats, stat)
		#two tailed
		if percentile >= 50:
			p = (1 - (percentile / 100)) * 2
		elif percentile < 50:
			p = (percentile / 100) * 2
		pvaluedict[index] = p
		condameandict[index] = condAmean
		condbmeandict[index] = condBmean
		condaMADdict[index] = condAmad
		condbMADdict[index] = condBmad

	#correct p values using BH
	positions = list(pvaluedict.keys())
	pvalues = list(pvaluedict.values())
	correctedpvalues = multipletests(pvalues, method = 'fdr_bh')[1]
	correctedpvalues = [float('{:.2e}'.format(fdr)) for fdr in correctedpvalues]
	correctedpvaluedict = dict(zip(positions, correctedpvalues))

	#turn into df
	pdf = pd.DataFrame.from_dict(correctedpvaluedict, orient = 'index')
	pdf.columns = ['FDR']

	condameandf = pd.DataFrame.from_dict(condameandict, orient = 'index')
	condameandf.columns = ['{0}_mean'.format(conditionA)]

	condbmeandf = pd.DataFrame.from_dict(condbmeandict, orient = 'index')
	condbmeandf.columns = ['{0}_mean'.format(conditionB)]

	condamaddf = pd.DataFrame.from_dict(condaMADdict, orient = 'index')
	condamaddf.columns = ['{0}_MAD'.format(conditionA)]

	condbmaddf = pd.DataFrame.from_dict(condbMADdict, orient = 'index')
	condbmaddf.columns = ['{0}_MAD'.format(conditionB)]

	bigdf = pd.concat([pdf, condameandf, condbmeandf, condamaddf, condbmaddf], axis = 1, join = 'inner')
	
	return bigdf

def getpvalues_LME(countdf, sampconds, conditionA, conditionB):
	#Use a linear mixed effects model to interrogate differences
	#between conditions

	#in countdf, indicies (rownames) are utr_nt and columns are samples
	#log2 transform count values

	pvaluedict = OrderedDict() #{utr_nt : pvalue}
	condameans = []
	condbmeans = []
	condamads = []
	condbmads = []
	condasds = []
	condbsds = []
	
	sampconddf = pd.read_csv(sampconds, sep = '\t', index_col = False, header = 0)
	colnames = list(sampconddf.columns)
	covariate_columns = [c for c in colnames if 'covariate' in c]

	condasamps = sampconddf.loc[sampconddf['condition'] == conditionA, 'sample'].tolist()
	condbsamps = sampconddf.loc[sampconddf['condition'] == conditionB, 'sample'].tolist()
	condasamps = [os.path.basename(samp) for samp in condasamps]
	condbsamps = [os.path.basename(samp) for samp in condbsamps]

	allsamps = condasamps + condbsamps
	samps = OrderedDict({'conditionA' : condasamps, 'conditionB' : condbsamps})

	print('Condition A samples: ' + (', ').join(condasamps))
	print('Condition B samples: ' + (', ').join(condbsamps))

	#Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
	currentutr = ''
	for index, row in countdf.iterrows():
		utr = index.split('_')[0]
		nt = int(index.split('_')[1])
		if utr != 'Gdf11':
			continue
		if utr != currentutr:
			print('Calculating pvalue for {0}...'.format(utr))
			currentutr = utr

		d = {}
		d['nt'] = [index] * len(allsamps)
		d['variable'] = allsamps
		
		values = [] #counts
		for cond in samps:
			for sample in samps[cond]:
				value = row[sample]
				values.append(value)

		#log2transform
		pc = 0.01
		values = [log2(v + pc) for v in values]

		d['value'] = values

		if covariate_columns:
			for covcol in covariate_columns:
				covs= [] #classifications for this covariate 
				for samp in allsamps:
					cov = sampconddf.loc[sampconddf['sample'] == samp, covcol].tolist()[0]
					covs.append(cov)
				d[covcol] = covs

			covariatestring = '+'.join(covariate_columns)

		condavalues = values[:len(condasamps)]
		condbvalues = values[len(condasamps):]
		condamean = np.mean(condavalues)
		condbmean = np.mean(condbvalues)
		condamad = MAD(condavalues, c = 1)
		condbmad = MAD(condbvalues, c = 1)
		condasd = np.std(condavalues)
		condbsd = np.std(condbvalues)

		condameans.append(condamean)
		condbmeans.append(condbmean)
		condamads.append(condamad)
		condbmads.append(condbmad)
		condasds.append(condasd)
		condbsds.append(condbsd)

		#If there is an NA count value, we are not going to calculate a pvalue for this gene
		p = None
		if np.isnan(values).any():
			p = np.nan

		conds = []
		for cond in samps:
			conds += [cond] * len(samps[cond])
		cond1s = []
		cond2s = []
		for cond in conds:
			if cond == 'conditionA':
				cond1s.append(1)
				cond2s.append(0)
			elif cond == 'conditionB':
				cond1s.append(0)
				cond2s.append(1)

		d['cond1'] = cond1s #e.g. [1, 1, 1, 1, 0, 0, 0, 0]
		d['cond2'] = cond2s #e.g. [0, 0, 0, 0, 1, 1, 1, 1]
		d['samples'] = [x + 1 for x in range(len(allsamps))]

		#Turn this dictionary into a DataFrame
		rowdf = pd.DataFrame.from_dict(d)

		#Get LME pvalue, but only if we haven't already determined that the pvalue is NA because we are missing one or more psi values
		#Lots of warnings about convergence, etc. Suppress them.
		if not p:
			with warnings.catch_warnings():
				warnings.filterwarnings('ignore')

				#So apparently, some combinations of psi values will give nan p values due to a LinAlgError that arises from a singular
				#hessian matrix during the fit of the model.  However, the equivalent code in R (nlme::lme) never gives this error, even with
				#the same data. It's not clear from just looking at the psi values why this is.  However, I found that by varying the 
				#start_params in the fit, this can be avoided. If this is done, the resulting p value always matches what is given in R.
				#Further, the p value is the same regardless of the start_param.
				#But it's not clear to me why changing the start_param matters, or what the default is here or with nlme.
				#So let's try a few starting paramters.  Regardless, this seems to affect a small number of genes (<1%), and it is causing 
				#false negatives because genes that should get p values (may or may not be sig) are getting NA.
				possible_start_params = [0, 1, -1, 2, -2, -3, 3]
				numberoftries = -1
				for param in possible_start_params:
					#if we already have a pvalue, don't try again
					if p != None and not np.isnan(p):
						break
					#First time through, numberoftries = 0, and we are just using a placeholder startparam (0) here because we aren't even using it.
					#Gonna use whatever the default is
					numberoftries +=1
					try:
						#actual model
						if covariate_columns:
							md = smf.mixedlm('value ~ cond1' + '+' + covariatestring, data = rowdf, groups = 'samples', missing = 'drop')
						else:
							md = smf.mixedlm('value ~ cond1', data = rowdf, groups = 'samples', missing = 'drop')
						if numberoftries == 0:
							mdf = md.fit(reml = False) #REML needs to be false in order to use log-likelihood for pvalue calculation
						elif numberoftries > 0:
							mdf = md.fit(reml = False, start_params = [param])

						#null model
						if covariate_columns:
							nullmd = smf.mixedlm('value ~ ' + covariatestring, data = rowdf, groups = 'samples', missing = 'drop')
						else:
							nullmd = smf.mixedlm('value ~ 1', data = rowdf, groups = 'samples', missing = 'drop')
						if numberoftries == 0:
							nullmdf = nullmd.fit(reml = False)
						elif numberoftries > 0:
							nullmdf = nullmd.fit(reml = False, start_params = [param])

						#Likelihood ratio
						LR = 2 * (mdf.llf - nullmdf.llf)
						p = chi2.sf(LR, df = 1)

					#These exceptions are needed to catch cases where either all psi values are nan (valueerror) or all psi values for one condition are nan (linalgerror)
					except (ValueError, np.linalg.LinAlgError):
						p = np.nan

		pvaluedict[index] = float('{:.2e}'.format(p))

	#Correct pvalues using BH method, but only using pvalues that are not NA
	pvalues = list(pvaluedict.values())
	pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'nan']
	fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
	fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

	#Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
	fdrswithnas = []
	fdrindex = 0
	for pvalue in pvalues:
		#print(pvalue)
		if str(pvalue) != 'nan':
			fdrswithnas.append(fdrs[fdrindex])
			fdrindex +=1
		elif str(pvalue) == 'nan':
			fdrswithnas.append(np.nan)


	fdrdict = dict(zip(list(pvaluedict.keys()), fdrswithnas))

	#Turn into a df
	pdf = pd.DataFrame.from_dict(fdrdict, orient = 'index')
	pdf.columns = ['FDR']
	#Add means and MADs
	pdf = pdf.assign(condamean = condameans)
	pdf = pdf.assign(condbmean = condbmeans)
	pdf = pdf.assign(condamad = condamads)
	pdf = pdf.assign(condbmad = condbmads)
	pdf = pdf.assign(condasd = condasds)
	pdf = pdf.assign(condbsd = condbsds)

	#rename columns
	pdf.columns = ['FDR', '{0}_mean'.format(conditionA), '{0}_mean'.format(conditionB), '{0}_MAD'.format(conditionA), '{0}_MAD'.format(conditionB), '{0}_sd'.format(conditionA), '{0}_sd'.format(conditionB)]

	return pdf


#TODO
#To define "regions":
#start at sig nt, allow max gap of <gap> nonsig nt

#SD of log2FC


#Get relationship of oligo coordinates to genome coordinates
oligo2genomedf = oligos2genome(sys.argv[4])

#Get df of umi counts
umidf = collatereplicates(sys.argv[1], sys.argv[2], sys.argv[3])

#Quantile normalize counts
umidf_qn = normalizecounts(umidf)
umidf.to_csv(path_or_buf = 'umicounts.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')
umidf_qn.to_csv(path_or_buf = 'umicounts_qn.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')

#Get nt level coverage as the median count across oligos that cover each nt
ntdf = getntcoverage(umidf_qn, 4, 260)
ntdf.to_csv(path_or_buf = 'umicounts_qn_nt.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')

#Smooth coverage with 11 nt window (5 nt +/- center of window)
smoothntdf = smoothcounts(ntdf, 11)
smoothntdf.to_csv(path_or_buf = 'umicounts_qn_nt_smooth.txt', sep = '\t', header = True, index = True, index_label = 'position', float_format = '%.3f')

smoothntdf = pd.read_csv('umicounts_qn_nt_smooth.txt', sep = '\t', header = 0, index_col = 'position')
#Calculate mean counts across replicates and log2FC across conditions
log2fcdf = getlog2fc(smoothntdf, sys.argv[1], sys.argv[2], sys.argv[3])

#Calculate p values at each nt using either MAD or LME
#pvaldf = getpvalues_MADstat(smoothntdf, sys.argv[1], sys.argv[2], sys.argv[3])
pvaldf = getpvalues_LME(smoothntdf, sys.argv[1], sys.argv[2], sys.argv[3])

#Merge log2FC df and pvalue df on position (utr_nt)
log2fc_and_pvaldf = pd.merge(log2fcdf, pvaldf, left_index = True, right_index = True)
log2fc_and_pvaldf['position'] = log2fc_and_pvaldf.index

#Split position column into 'Gene' and 'utrpos'
s = log2fc_and_pvaldf['position'].str.split('_', n = 1)
log2fc_and_pvaldf['Gene'], log2fc_and_pvaldf['utrcoord'] = s.str[0], s.str[1]
cols = log2fc_and_pvaldf.columns.tolist()
colnames = [cols[-2], cols[-1]]
colnames += cols[0:9]
log2fc_and_pvaldf = log2fc_and_pvaldf[colnames]
#Set utrcoord to int64 (was object)
log2fc_and_pvaldf['utrcoord'] = log2fc_and_pvaldf.utrcoord.astype(int)
log2fc_and_pvaldf = log2fc_and_pvaldf.drop(['position'], axis = 1)

#Merge genomecoord df and log2fc_and_pvaldf
newdf = pd.merge(oligo2genomedf, log2fc_and_pvaldf, how = 'inner', on = ['Gene', 'utrcoord'])
print(newdf.head())
newdf.to_csv(path_or_buf = 'mpraresults.txt', sep = '\t', header = True, index = False, float_format = '%.3g')
