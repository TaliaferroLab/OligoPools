#Given a gff annotation of a genome and a genome sequence, make a pool of oligo sequences that cover all 3' UTRs.
#This can be subsetted to only consider some genes (e.g. only genes expressed in certain cells).
#Also, the parameters for oligo coverage (oligo length, coverage depth, etc.) can be tuned.
#For a start, each oligo is a 110mer and the average coverage is 2.5X

import gffutils
import os
import sys
from operator import itemgetter
from itertools import groupby
import subprocess
import pybedtools
import pickle
import pandas as pd
from Bio import SeqIO
import gzip
import argparse
from biomaRt import getmouseens2dmelens
from rpy2.rinterface import RRuntimeError
import collections

#Write UTR coords to file so that we can use Bedtools.
#geneUTRcoords is a list of lists.
def writeUTRcoords(geneUTRcoords):
	with open('temp.bed', 'w') as outfh:
		for exon in geneUTRcoords:
			outfh.write(('\t').join(exon) + '\n')

	#Now sort the bed by chrm and start position
	command = ['sort', '-k1,1', '-k2,2n', 'temp.bed']
	with open('temp.sorted.bed', 'w') as outfh:
		subprocess.call(command, stdout = outfh)


#Get a set of polyAsite coords from another genome.  This will be used as a reference for determining whether a given polyAsite
#is conserved.  Write these to file, then use UCSC liftover to convert them to another genome's coordinate.
#In the end, return a dictionary of "valid" positions for polyAsites in the liftedover genome.  In this example, the supplied gff would be
#for hg38, and the returned coordinates of valid positions would be in mm10.

def getothergenomepolyA(gff):
	#Not looking for orthologous polyA in drosophila
	return None

	#Make gff database

	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')

	polyAsites = {} #{transcriptid : [chrm, start, stop, strand]}
	genecounter = 0

	for gene in genes:
		if gene.attributes['gene_type'][0] != 'protein_coding':
			continue
		genecounter +=1
		if genecounter % 5000 == 0:
			print 'Gene {0}...'.format(genecounter)

		for transcript in db.children(gene, featuretype = 'transcript'):
			#Only look at protein coding transcripts
			if transcript.attributes['transcript_type'][0] != 'protein_coding':
				continue
			tags = transcript.attributes['tag']
			#Don't consider transcripts with these tags
			if 'cds_start_NF' in tags or 'mRNA_start_NF' in tags or 'cds_end_NF' in tags or 'mRNA_end_NF' in tags:
				continue

			strand = transcript.strand

			if strand == '+':
				polyAsite = transcript.end
			elif strand == '-':
				polyAsite = transcript.start
			polyAsites[str(transcript.id)] = [str(transcript.chrom), polyAsite, polyAsite + 1, transcript.strand]

	#Write the polyAsites to a bed file.
	with open('polyAsites.bed', 'w') as outfh:
		for transcript in polyAsites:
			chrm = polyAsites[transcript][0]
			start = polyAsites[transcript][1]
			stop = polyAsites[transcript][2]
			strand = polyAsites[transcript][3]
			outfh.write(('\t').join([chrm, str(start), str(stop), strand, transcript]) + '\n')

	#Liftover these polyA sites to mm10
	command = ['liftOver', 'polyAsites.bed', '/Users/mtaliaferro/Desktop/Annotations/hg38ToMm10.over.chain.gz', 'liftedpolyAsites.bed', 'unliftedpolyAsites.bed']
	subprocess.call(command)
	num_lifted = sum(1 for line in open('liftedpolyAsites.bed'))
	num_unlifted = sum(1 for line in open('unliftedpolyAsites.bed')) / 2 #In unliftedpolyAsites.bed, there are 2 lines per polyAsite
	percentlifted = round((num_lifted / float(num_lifted + num_unlifted)) * 100, 2)
	print 'Lifted over {0} percent of {1} polyAsites.'.format(percentlifted, num_unlifted + num_lifted)
	os.remove('polyAsites.bed')

	#Define "valid" polyAsite regions.
	#Within 200 nt of a lifted over site.  Strand does not matter.
	validwindows = {} #{chrm : [list of valid nt]}
	with open('liftedpolyAsites.bed', 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			chrm = line[0]
			start = int(line[1])
			windowstart = start - 200 #Can change this
			windowstop = start + 200 #Can change this
			if chrm not in validwindows:
				validwindows[chrm] = []
			windowrange = range(windowstart, windowstop + 1)
			validwindows[chrm] += windowrange

	#Turn list of valid positions into sets for faster searching
	for chrm in validwindows:
		validwindows[chrm] = set(validwindows[chrm])

	os.remove('liftedpolyAsites.bed')
	os.remove('unliftedpolyAsites.bed')
	return validwindows

def getgenenames(gff):
	#Go through the gff and get gene short names
	print 'Getting gene names...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)

	genes = db.features_of_type('gene')

	ens2genename = {} #{ensid : gene short name}

	for gene in genes:
		genename = str(gene.id).split(':')[1]
		shortname = gene.attributes['Name'][0]
		ens2genename[genename] = shortname

	print 'Found names for {0} genes.'.format(len(ens2genename))

	return ens2genename

#Iterate through a gff, gene by gene, collecting all 3' UTR coords for protein coding transcripts of that gene
def iterategff(gff, validwindows):
	#Make gff database

	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	genes = db.features_of_type('gene')
	genecounter = 0
	geneswithoutpassingtranscripts = []
	transcriptcounter = 0
	longUTRcounter = 0
	tslcounter = 0
	notinvalidwindowscounter = 0

	mergedexons = {} #{genename : [chrm, strand, [[mergedexon1start, mergedexon1stop], [mergedexon2start, mergedexon2stop]]]}

	for gene in genes:
		genename = str(gene.id).split(':')[1]

		chrm = str(gene.chrom)
		if gene.attributes['biotype'][0] != 'protein_coding':
			continue
		genecounter +=1
		if genecounter % 5000 == 0:
			print 'Gene {0}...'.format(genecounter)

		#Does any UTR in this gene have introns
		hasjunctions = False


		#First get all coordinates that are UTR coordinate for any transcript of this gene.
		geneUTRcoords = []
		for transcript in db.children(gene, featuretype = 'transcript'):
			
			#Only look at protein coding transcripts
			if transcript.attributes['biotype'][0] != 'protein_coding':
				continue

			#Drosophila annotation doesn't have tsl or tags
			'''
			try:
				tsl = transcript.attributes['transcript_support_level'][0]
			except KeyError: #if the transcript doesn't have a transcript_support_level annotation
				continue
			
			tags = transcript.attributes['tag']
			#Don't consider transcripts with these tags
			if 'cds_start_NF' in tags or 'mRNA_start_NF' in tags or 'cds_end_NF' in tags or 'mRNA_end_NF' in tags:
				continue
			'''

			transcriptcounter +=1

			#Only consider transcripts with transcript_support_level of 1
			#Drosophila annotation doesn't have tsl
			'''
			if tsl == 'NA' or int(tsl) > 1:
				tslcounter +=1
				#print 'TSL fail, {0}.'.format(str(transcript.id))
				continue
			'''

			#If the polyA site for this transcript is not a "valid" polyA site as defined by validwindows, forget it
			strand = transcript.strand
			if strand == '+':
				polyAsite = transcript.end
			elif strand == '-':
				polyAsite = transcript.start

			#Not looking for orthologous polyA for drosophila
			'''
			if polyAsite not in validwindows[chrm]:
				notinvalidwindowscounter +=1
				continue
			'''

			exoncoords = [] #[list of coords that are exonic]
			CDScoords = [] #[list of coords that are CDS exonic]

			for exon in db.children(transcript, featuretype = 'exon', order_by = 'start'):
				exoncoords += range(exon.start, exon.end + 1)
			for CDSexon in db.children(transcript, featuretype = 'CDS', order_by = 'start'):
				CDScoords += range(CDSexon.start, CDSexon.end + 1)

			#3' UTR starts directly after the CDS end
			if strand == '+':
				UTRstart = max(CDScoords) + 1
				UTRend = transcript.end
			elif strand == '-':
				UTRstart = transcript.start
				UTRend = min(CDScoords) - 1

			#Get 64 nt immediately upstream of the beginning of the UTR
			#If the UTR is shorter than 96 (160 - one step size), we need extra sequence to fit one whole oligo.
			#The amount of extra sequence we need 160 - 64 - UTRlength

			#Get 44 nt immediately  upstream of the beginning of the UTR.
			#If the UTR is shorter than 66 (110 - one step size), we need extra sequence to fit one whole oligo.
			#The amount of extra sequence we need is 110 - 44 - UTRlength
			if strand == '+':
				UTRstartindex = exoncoords.index(max(CDScoords))
				#upstreamnt = exoncoords[UTRstartindex - 44 : UTRstartindex + 1]
				upstreamnt = exoncoords[UTRstartindex - 43 : UTRstartindex + 1]
				#upstreamnt = exoncoords[UTRstartindex - 63 : UTRstartindex + 1]
				UTRl = (UTRend - UTRstart) + 1
				if UTRl < 66:
					extrantneeded = 66 - UTRl
					upstreamnt = exoncoords[UTRstartindex - 43 - extrantneeded : UTRstartindex + 1]
					#print UTRl, len(upstreamnt), strand 
			elif strand == '-':
				UTRstartindex = exoncoords.index(min(CDScoords))
				#upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 45]
				upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 44]
				#upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 64]
				UTRl = (UTRend - UTRstart) + 1
				if UTRl < 66:
					extrantneeded = 66 - (UTRend - UTRstart + 1)
					upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 44 + extrantneeded]
					#print UTRl, len(upstreamnt), strand 

			UTRcoords = range(UTRstart, UTRend + 1)
			#Add on the upstream 64 nt
			UTRcoords = upstreamnt + UTRcoords

			#Filter for those that are exonic
			UTRcoords = sorted(list(set(UTRcoords).intersection(exoncoords)))

			#If this UTR is longer than 8 kb or shorter than 96 (160 - one step size) nt, forget it
			if len(UTRcoords) > 100000000:
				#print 'Long UTR : {0}'.format(str(transcript.id))
				longUTRcounter +=1
				continue

			#Now get breaks in consecutive exonic positions
			#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
			UTRexoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
			for k, g in groupby(enumerate(UTRcoords), lambda (index, item): index-item):
				exonbp = map(itemgetter(1), g)
				if len(exonbp) > 1:
					UTRexoncoords.append([exonbp[0], exonbp[-1]])

			#Add this exon(s) to geneUTRcoords
			for exon in UTRexoncoords:
				geneUTRcoords.append([gene.chrom, str(exon[0]), str(exon[1])])

			#If this UTR has more than one exon, we need to make junction oligos for every exon/exon boundary
			if len(UTRexoncoords) == 1:
				pass
			elif len(UTRexoncoords) > 1:
				hasjunctions = True
				junctionblocks = [] #This is the entire region covering a junction, from 99 nt upstream to 99 nt downstream
				for i in range(len(UTRexoncoords) - 1):
					lastexonicposition = UTRexoncoords[i][1] #last exonic coord (left to right) before the intron
					firstexonicposition = UTRexoncoords[i + 1][0] #first exonic coord (left to right) after the intron
					#If this junction point is upstream of where the UTR starts (the real UTR start, not the part with the one step size added on), then there's
					#no reason to make a junction oligo for this junction.  Basically, I don't want to be making junction oligos for CDS junctions.
					if strand == '+':
						if firstexonicposition < UTRstart:
							continue
					elif strand == '-':
						if firstexonicposition > UTRend:
							continue
					lastexonicindex = exoncoords.index(lastexonicposition)
					firstexonicindex = exoncoords.index(firstexonicposition)
					#Assuming an oligo length of 110, the evenly spaced bridge should be 54upstream--last exonic--first exonic--54 downstream
					#Upstream tiled should be offset by 44 from above: 98upstream--lastexonic--firstexonic--10downstream
					#Downstream tiled should be offset by 44 the other way: 10upstream--lastexonic--firstexonic--98downstream
					#So what we really want is 98upstream--lastexonic--firstexonic--98downstream as a block.  Then make oligos from that later.
					#If 98 nt away from the junction is outside the bounds of the transcript, add sequence from beyond the transcript to make up the difference.

					#Assuming an oligo length of 160, the evenly spaced bridge should be 79upstream--last exonic--first exonic--79 downstream
					#Upstream tiled should be offset by 64 from above: 143upstream--lastexonic--firstexonic--15downstream
					#Downstream tiled should be offset by 64 the other way: 15upstream--lastexonic--firstexonic--143downstream
					#So what we really want is 143upstream--lastexonic--firstexonic--143downstream as a block.  Then make oligos from that later.
					#If 98 nt away from the junction is outside the bounds of the transcript, add sequence from beyond the transcript to make up the difference.

					'''
					if lastexonicindex - 143 < 0: #If the junction is within 143 nt of the transcript start, 
						missingntleft = 144 - lastexonicindex
						junctionblockstart = exoncoords[0] - missingntleft
						junctionblock = range(junctionblockstart, exoncoords[0]) + exoncoords[0:lastexonicindex + 1] + exoncoords[firstexonicindex : firstexonicindex + 143 + 1]

					elif firstexonicindex + 143 > len(exoncoords) - 1: #if the junction is within 143 nt of the transcript end
						rightntintranscript = len(exoncoords[firstexonicindex:])
						missingntright = 144 - rightntintranscript ###CHECK
						junctionblockend = exoncoords[-1] + missingntright
						junctionblock = exoncoords[lastexonicindex - 143 : lastexonicindex + 1] + exoncoords[firstexonicindex :] + range(exoncoords[-1] + 1, junctionblockend + 1)

					elif lastexonicindex - 143 >=0 and firstexonicindex + 143 <= len(exoncoords) - 1: #if the entire junction block is contained within the transcript
						junctionblockstart = lastexonicindex - 143
						junctionblockend = firstexonicindex + 143 + 1
						junctionblock = exoncoords[junctionblockstart : junctionblockend]
					'''

					if lastexonicindex - 98 < 0: #If the junction is within 98 nt of the transcript start, 
						missingntleft = 99 - lastexonicindex
						junctionblockstart = exoncoords[0] - missingntleft
						junctionblock = range(junctionblockstart, exoncoords[0]) + exoncoords[0:lastexonicindex + 1] + exoncoords[firstexonicindex : firstexonicindex + 98 + 1]

					elif firstexonicindex + 98 > len(exoncoords) - 1: #if the junction is within 143 nt of the transcript end
						rightntintranscript = len(exoncoords[firstexonicindex:])
						missingntright = 99 - rightntintranscript
						junctionblockend = exoncoords[-1] + missingntright
						junctionblock = exoncoords[lastexonicindex - 98 : lastexonicindex + 1] + exoncoords[firstexonicindex :] + range(exoncoords[-1] + 1, junctionblockend + 1)

					elif lastexonicindex - 98 >=0 and firstexonicindex + 98 <= len(exoncoords) - 1: #if the entire junction block is contained within the transcript
						junctionblockstart = lastexonicindex - 98
						junctionblockend = firstexonicindex + 98 + 1
						junctionblock = exoncoords[junctionblockstart : junctionblockend]					

					junctionblocks.append(junctionblock)


		#Write geneUTRcoords to file so we can use Bedtools
		writeUTRcoords(geneUTRcoords)

		#Merge the UTR exons
		bed = pybedtools.BedTool('temp.sorted.bed')
		mergedbed = bed.merge()

		#Get the merged "exons"
		exons = []
		for exon in mergedbed:
			exons.append([int(exon[1]), int(exon[2])])

		if not exons:
			#It is possible that this gene has protein-coding transcripts but that none of them pass all the filters
			#(TSL > 1, polyA site within 200 nt of human polyA site, no tag like 'mRNA_end_NF', UTR less than 8kb)
			#If this is the case, then don't put this gene (and it's empty exons) into mergedexons
			geneswithoutpassingtranscripts.append(genename)
			continue

		#Put the merged "exons" and the junction blocks in mergedexons
		if hasjunctions == True:
			mergedexons[genename] = [gene.chrom, gene.strand, exons, junctionblocks]
		elif hasjunctions == False:
			mergedexons[genename] = [gene.chrom, gene.strand, exons, []]

	print 'Looked through {0} protein coding genes. {1} of these ({2}%) had at least one transcript that passed all filters.'.format(genecounter, genecounter - len(geneswithoutpassingtranscripts), round(((genecounter - len(geneswithoutpassingtranscripts)) / float(genecounter)) * 100), 3)
	print 'Looked through {0} transcripts. {1} had UTRs that were too long or too short. {2} were not supported in the annotation (TSL > 1). {3} did not have polyA sites within 200 nt of the orthologous human genome location.'.format(transcriptcounter, longUTRcounter, tslcounter, notinvalidwindowscounter)
	os.remove('temp.bed')
	os.remove('temp.sorted.bed')
	
	return mergedexons


def filterexpression(tpmtable, mergedexons):
	print 'Filtering by expression...'
	genesbeforefilter = len(mergedexons)
	#Given a table of tpm expression values, filter for genes that pass some expression cutoff.
	df = pd.read_csv(tpmtable, sep = '\t', header = 0)

	'''
	#The P70brain1 sample looked funky by PCA/clustering, so I dropped it
	#Also, the cerebellum is quite different than the rest of the brain, so I dropped those samples too

	df = df.drop(['P70brain1', 'P56cerebellum1', 'P56cerebellum2'], axis = 1)

	#Add a column of the minimum expression across all samples for that gene
	df['minexp'] = df.loc[:, ['E18brain1', 'E18brain2', 'P56cortex1', 'P56cortex2', 'P56cortex3', 'P56cortex4', 'P70brain2']].min(axis = 1)
	'''
	
	#Get all genes that are expressed at least x tpm in ANY sample
	#df = df.query('maxexp >= 10')
	#Get all genes that are expressed at least 20 tpm in all head samples
	df = df.query('passexp == \'yes\'')
	expressedgenes = df['Gene'].tolist()

	#For genes in mergedexons, if that gene isn't in expressedgenes, remove it
	for gene in mergedexons.keys():
		if gene not in expressedgenes:
			del mergedexons[gene]

	print '{0} of {1} genes pass expression filter.'.format(len(mergedexons), genesbeforefilter)
	return mergedexons

def filterlocalization(lrtable, mergedexons):
	#In this strategy, instead of filtering for genes that are expressed in brain, filter for genes that are known to be localized.
	#This way we will find the fragment within the 3' UTR that is driving localization.
	#The LR data comes from the "Bartel analysis" and is usually FractionationLR_z.txt.
	#It's z scores of LR values for almost every fractionation/sequencing experiment we've done (about 30 of them).

	#For the human oligos, this LR is still LR of mouse genes.  What we need to get is the human orthologs of the 
	#mouse genes that pass the LR filter.

	print 'Filtering by localization...'
	genesbeforefilter = len(mergedexons)
	#Given a table of tpm expression values, filter for genes that pass some expression cutoff.
	df = pd.read_csv(lrtable, sep = '\t', header = 0)

	df = pd.read_csv(lrtable, sep = '\t', header = 0)
	#Get rid of the 'ribo' and 'mito' columns
	#Also get rid of the Mbnl2_1 experiment. It's an outlier to every other experiment here.
	df = df.drop(['ribo', 'mito', 'Mbnl2_1'], axis = 1)
	#Get rid of any row that has an NA (didn't pass neurite and/or soma expression filter in a sample (tpm >= 5))
	df = df.dropna(axis = 0, how = 'any')
	#Get the median of all LRs for a gene
	med = df.median(axis = 1, skipna = True, numeric_only = True)
	df = df.assign(median = med)
	#Filter for those that have <= -0.75 LRz (negctrls), > 0.5 LRz, or are Ranbp1 or Actb
	query = 'median <= -0.7 | median > 0.38 | Gene == "Actb" | Gene == "Ranbp1"'
	df = df.query(query)
	chosengenes = df['ensembl_gene_id'].tolist()

	print '{0} mouse genes pass LR filters.'.format(len(chosengenes))

	#Get mouse/human ortholog relationships
	#biomaRt is a bitch and won't connect all the time anymore
	#try to connect until we do
	hasm2d = False
	while hasm2d == False:
		try:
			m2d = getmouseens2dmelens()
			hasm2d = True
		except RRuntimeError:
			print 'Trying to reconnect to biomaRt...' #try again

	#See which human orthologs pass expression
	dmelpassinggenes = []
	for gene in chosengenes:
		try:
			dmelgene = m2d[gene]
		except KeyError:
			#print 'Mouse gene {0} not in mouse2dmel dictionary.'.format(gene)
			pass
		dmelpassinggenes.append(dmelgene)

	dmelpassinggenes = list(set(dmelpassinggenes))
	print 'Starting with {0} mouse genes that pass LR filter, found {1} Drosophila orthologs.'.format(len(chosengenes), len(dmelpassinggenes))

	#For genes in mergedexons, if that gene isn't in humanpassinggenes, remove it
	for gene in mergedexons.keys():
		if gene not in dmelpassinggenes:
			del mergedexons[gene]

	print '{0} of {1} genes pass LR filter.'.format(len(mergedexons), genesbeforefilter)
	return mergedexons

def breakoligo(coords):
	#Get breaks in consecutive oligo positions
	#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
	brokencoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
	for k, g in groupby(enumerate(coords), lambda (index, item): index-item):
		exonbp = map(itemgetter(1), g)
		if len(exonbp) > 1:
			brokencoords.append([exonbp[0], exonbp[-1]])

	return brokencoords

def cutexons(coords):
	#For an entry in mergedexons, split it up into oligos
	strand = coords[1]
	exons = coords[2]
	junctions = coords[3]
	oligos = [] #[[oligo1start, oligo1end], [oligo2start, oligo2end]]
	junctionoligos = [] #[[j1half1start, j1half1end], [j1half2start, j1half2end], [j2half1start, j2half1end], [j2half2start, j2half2end]]

	#If there's only one exon, this is easy.
	if len(exons) == 1:
		if strand == '+':
			for exon in exons:
				exonstart = exon[0]
				exonend = exon[1]
				currentpos = exonstart
				#This is confusing but because the gff is 1-based and open, an oligo that is x to x + 109 will actually be length 110
				#This is confusing but because the gff is 1-based and open, an oligo that is x to x + 159 will actually be length 160
				while currentpos + 109 <= exonend:
					oligostart = currentpos
					oligoend = currentpos + 109
					oligos.append([[oligostart, oligoend]])
					currentpos += 44 #for 2.5X coverage (160 / 2.5)
				#If the final oligo does not lie flush with the end of the exon, make one more oligo that does
				if currentpos != exonend:
					oligos.append([[exonend - 109, exonend]])

		elif strand == '-':
			#Gotta start from the other end of coords so that we make the same kind of "final" oligo that lies flush with the polyA site
			#If we didn't consider strand and treated the two strands the same, the final flush oligo for the - strand would actually be flush with the stop codon
			#or 50 nt upstream of the stop codon where the exon starts
			exons = reversed(exons)
			for exon in exons:
				exonstart = exon[1] #these coords have start > end because they are on minus strand
				exonend = exon[0]
				currentpos = exonstart
				while currentpos - 109 >= exonend:
					oligoend = currentpos #these coords are back to having start < end
					oligostart = currentpos - 109
					oligos.append([[oligostart, oligoend]])
					currentpos -= 44 #for 2.5X coverage (160 / 2.5)
				#If the final oligo does not lie flush with the end of the exon, make one more oligo that does
				if currentpos != exonend:
					oligos.append([[exonend, exonend + 109]])

	#If there's more than one exon, it's more complicated
	elif len(exons) > 1:
		if strand == '+':
			exoniccoords = []
			for exon in exons:
				exoniccoords += range(exon[0], exon[1] + 1)
			currentposindex = 0
			#This is slightly different than the single-exon example above. Here we are stepping along indices, so we need the oligo to go from index x to index x + 160.
			#This will give coordinates that are 159 "places" apart
			while currentposindex + 110 <= len(exoniccoords):
				oligostart = currentposindex
				oligoend = currentposindex + 110
				oligocoords = exoniccoords[oligostart : oligoend]
				#We might need to break this oligo because it may have crossed an exon/exon boundary
				oligocoords = breakoligo(oligocoords) #if this doesn't cross a boundary, its now [[start, stop]]. if it does, its [[start1, stop1], [start2, stop2], [start3, stop3]]
				oligos.append(oligocoords) 
				currentposindex += 44
			#If the final oligo does not lie flush with the end of the last exon, make one more oligo that does
			if currentposindex != len(exoniccoords):
				oligocoords = exoniccoords[-110:]
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)

		elif strand == '-':
			exoniccoords = []
			for exon in exons:
				exoniccoords += range(exon[0], exon[1] + 1)
			exoniccoords.reverse()
			currentposindex = 0
			while currentposindex + 110 <= len(exoniccoords):
				oligostart = currentposindex + 110 #have to make it so start < end
				oligoend = currentposindex
				oligocoords = exoniccoords[oligoend : oligostart]
				#We might need to break this oligo because it may have crossed an exon/exon boundary
				oligocoords.reverse() #flip it back around so that startcoord < endcoord
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)
				currentposindex += 44
			#If the final oligo does not lie flush with the end of the last exon, make one more oligo that does
			if currentposindex != len(exoniccoords):
				oligocoords = exoniccoords[-110:]
				oligocoords = reversed(oligocoords)
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)


	for junction in junctions:
		jo1 = junction[0:110]
		jo2 = junction[44:154]
		jo3 = junction[88:198]

		#jo1 = junction[0:160]
		#jo2 = junction[64:224]
		#jo3 = junction[128:288]

		for oligo in [jo1, jo2, jo3]:
			brokencoords = breakoligo(oligo)
			junctionoligos.append(brokencoords)

	return oligos, junctionoligos

def makeoligos(mergedexons):
	#Go through merged exons.
	#Feed every gene to cutexons to receive back oligo coordinates
	oligocoords = {} #{gene : [chrm, strand, [oligos]]} #For oligos that do not have introns, i.e. continuous oligos
	junctionoligocoords = {} #{gene : [chrm, strand, [junction oligos]]} #For oligos that do have introns

	for gene in mergedexons:
		coords = mergedexons[gene]
		chrm = coords[0]
		strand = coords[1]
		oligos, junctionoligos = cutexons(coords)
		if oligos:
			oligocoords[gene] = [chrm, strand, oligos]
		if junctionoligos:
			junctionoligocoords[gene] = [chrm, strand, junctionoligos]

	return oligocoords, junctionoligocoords

def flatten(x):
	#Recursively flatten a list
	#https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists?page=1&tab=votes#tab-top
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]


def makegff(oligocoords, junctionoligocoords, genenames):
	with open('oligos.dm6.gff', 'w') as outfh:
		for gene in oligocoords:
			chrm = oligocoords[gene][0]
			strand = oligocoords[gene][1]
			oligos = oligocoords[gene][2]
			oligooutlists = []
			if strand == '+':
				#UTRstart = oligos[0][0][0] #this is actually 50 nt upstream of the UTR start, but it's the start of the oligo-producing region
				#UTRend = oligos[-1][-1][1]
				boundaries = flatten(oligos)
				UTRstart = min(boundaries)
				UTRend = max(boundaries)
			elif strand == '-':
				#UTRstart = oligos[-1][-1][0] #start has to be less than end for a gff
				#UTRend = oligos[0][0][1]
				boundaries = flatten(oligos)
				UTRstart = min(boundaries)
				UTRend = max(boundaries)
			oligocounter = 0

			#Make lines for "regular" (non-junction) oligos
			for oligo in oligos:
				if len(oligo) == 1: #if this is a single-exon oligo
					oligostart = oligo[0][0]
					oligoend = oligo[0][1]
					#Somehow a small number of oligos (~10) are ending up at length 159 (pre adapter addition). I can't figure out why.  Fix this lazily here.
					oligolength = (oligoend - oligostart) + 1
					if oligolength < 109:
						continue
					elif oligolength == 109:
						if strand == '+':
							oligoend += 1
						elif strand == '-':
							oligostart -= 1
					oligocounter +=1
					oligooutlist = [chrm, 'dm6', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=regular_oneexon' + ';Parent={0}.UTR'.format(gene)]
					oligooutlists.append(oligooutlist)
				elif len(oligo) > 1: #if this oligo spans multiple exons
					oligocounter +=1
					numberofpieces = len(oligo)
					oligostart = oligo[0][0]
					oligoend = oligo[-1][1]
					oligooutlist = [chrm, 'dm6', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=regular_multiexon' + ';Parent={0}.UTR'.format(gene)]
					piececounter = 0
					oligooutlists.append(oligooutlist)
					for piece in oligo:
						piececounter +=1
						piecestart = piece[0]
						pieceend = piece[1]
						oligooutlist = [chrm, 'dm6', 'oligopiece', str(piecestart), str(pieceend), '.', strand, '.', 'ID={0}.{1}.{2}'.format(gene, oligocounter, piececounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}.{2}'.format(gene, oligocounter, piececounter) + ';oligo_type=regularmultiexon' + ';Parent={0}.{1}'.format(gene, oligocounter)]
						oligooutlists.append(oligooutlist)

			#Make lines for junction oligos
			if gene in junctionoligocoords:
				junctionoligos = junctionoligocoords[gene][2]
				junctioncounter = 0
				for oligo in junctionoligos:
					oligostart = oligo[0][0]
					oligoend = oligo[-1][1]
					oligocounter +=1
					junctioncounter +=1
					piececounter = 0
					oligooutlist = [chrm, 'dm6', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=junction' + ';Parent={0}.UTR'.format(gene)]
					oligooutlists.append(oligooutlist)
					for piece in oligo:
						piececounter +=1
						piecestart = piece[0]
						pieceend = piece[1]
						oligooutlist = [chrm, 'dm6', 'junctionpiece', str(piecestart), str(pieceend), '.', strand, '.', 'ID={0}.{1}.{2}'.format(gene, oligocounter, piececounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=junction' + ';Parent={0}.{1}'.format(gene, oligocounter)]
						oligooutlists.append(oligooutlist)

			#Make gene-level line
			#This isn't actually the UTR.  It includes 50 nt or so of the coding sequence. It may not encompass the entire upstream coding region as a junction oligo piece may be further upstream than the first regular oligo start.
			geneline = [chrm, 'dm6', 'UTR', str(UTRstart), str(UTRend), '.', strand, '.', 'gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';number_of_oligos={0}'.format(oligocounter) + ';ID={0}.UTR'.format(gene)]

			outfh.write(('\t').join(geneline) + '\n')
			for oligooutlist in oligooutlists:
				outfh.write(('\t').join(oligooutlist) + '\n')

def gff2fasta(gff, genomefasta):
	print 'Indexing gff...'
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print 'Done indexing!'

	print 'Indexing genome sequence...'
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta), 'fasta'))
	print 'Done indexing!'

	oligos = db.features_of_type('oligo')
	seqs = {} #{oligoname : seq}
	oligonames = [] #list of oligoIDs, this is so we can order the oligos when we write the fasta
	adapter5 = 'GGCTACTACGCGTGTTTGAG'
	adapter3 = 'ACCTGCAGGCGAAGCTCTAC'

	for oligo in oligos:
		chrm = oligo.chrom
		oligostart = oligo.start
		oligoend = oligo.end
		strand = oligo.strand
		oligotype = oligo.attributes['oligo_type'][0]
		oligo_id = oligo.attributes['oligo_id'][0]
		genename = oligo.attributes['gene_name'][0]

		if oligotype == 'regular_oneexon':
			oligoname = oligo_id + '|' + genename
			if strand == '+':
				seq = seq_dict[chrm].seq[oligostart - 1 : oligoend].upper() #gff is one-based, seq_dict is 0-based
			elif strand == '-':
				seq = seq_dict[chrm].seq[oligostart - 1 : oligoend].reverse_complement().upper()

			seqs[oligoname] = str(seq)
			oligonames.append(oligoname)

		elif oligotype == 'regular_multiexon':
			oligoname = oligo_id + '|' + genename
			oligopieceseqs = []
			for oligopiece in db.children(oligo, featuretype = 'oligopiece'):
				if strand == '+':
					seq = seq_dict[chrm].seq[oligopiece.start - 1 : oligopiece.end].upper()
				elif strand == '-':
					seq = seq_dict[chrm].seq[oligopiece.start - 1 : oligopiece.end].reverse_complement().upper()
				oligopieceseqs.append(seq)

			oligopieceseqs = [str(oligopieceseq) for oligopieceseq in oligopieceseqs]
			if strand == '-':
				oligopieceseqs.reverse()
			seq = ''.join(oligopieceseqs)
			seqs[oligoname] = seq
			oligonames.append(oligoname)		

		elif oligotype == 'junction':
			oligoname = oligo_id + '|' + genename
			junctionpieceseqs = []
			for junctionpiece in db.children(oligo, featuretype = 'junctionpiece'):
				if strand == '+':
					seq = seq_dict[chrm].seq[junctionpiece.start - 1 : junctionpiece.end].upper()
				elif strand == '-':
					seq = seq_dict[chrm].seq[junctionpiece.start - 1 : junctionpiece.end].reverse_complement().upper()
				junctionpieceseqs.append(seq)

			junctionpieceseqs = [str(junctionpieceseq) for junctionpieceseq in junctionpieceseqs]
			if strand == '-':
				junctionpieceseqs.reverse()
			seq = ''.join(junctionpieceseqs)
			seqs[oligoname] = seq
			oligonames.append(oligoname)

	with open('oligos.dm6.fa', 'w') as outfh:
		for oligoname in oligonames:
			seq = seqs[oligoname]
			seq = adapter5 + seq + adapter3
			outfh.write('>' + oligoname + '\n' + seq + '\n')

	os.remove(db_fn)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF file of genome to make oligos from.')
	#parser.add_argument('--othergenomegff', type = str, help = 'GFF file of genome in which we will look for orthologous polyA sites.')
	parser.add_argument('--filter', type = str, choices = ['expression', 'localization'], help = 'How to pick the genes to make oligos from? Brain-expressed genes or localized genes?')
	parser.add_argument('--filterdata', type = str, help = 'Data to perform the filtering with.  TPMtable.txt if expression. FractionationLR_z if localization.')
	parser.add_argument('--genomefasta', type = str, help = 'Sequence in fasta format for the genome to make oligos from.')
	args = parser.parse_args()

	#Find orthologous polyA sites
	#validwindows = getothergenomepolyA(args.othergenomegff)
	
	#Find UTRs that pass transcript and polyA site filters
	#mergedexons = iterategff(args.gff, validwindows)
	mergedexons = iterategff(args.gff, None)
	
	#Filter for expressed or localized genes
	if args.filter == 'expression':
		mergedexons = filterexpression(args.filterdata, mergedexons)
	elif args.filter == 'localization':
		mregedexons = filterlocalization(args.filterdata, mergedexons)
	#Get ensembl gene ID / genename relationships
	genenames = getgenenames(args.gff)


	#Print the total amount of UTR space we end up with
	UTRsize = 0
	for gene in mergedexons:
		for exon in mergedexons[gene][2]:
			UTRsize += exon[1] - exon[0]
		junctionblocks = mergedexons[gene][3]
		if junctionblocks:
			for junctionblock in junctionblocks:
				UTRsize += len(junctionblock)

	print 'After filtering, end up with {0} nt of UTR sequence.'.format(UTRsize)

	oligocoords, junctionoligocoords = makeoligos(mergedexons)

	#Print how many oligos this gets broken into
	oligocount = 0
	for gene in oligocoords:
		oligos = oligocoords[gene][2]
		oligocount += len(oligos)
	for gene in junctionoligocoords:
		oligos = junctionoligocoords[gene][2]
		oligocount += len(oligos)

	print 'This is broken up into {0} oligos.'.format(oligocount)

	#Make oligo gff
	makegff(oligocoords, junctionoligocoords, genenames)
	#Make fasta
	gff2fasta('oligos.dm6.gff', args.genomefasta)

