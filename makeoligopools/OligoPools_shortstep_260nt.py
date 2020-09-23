#Given a gff annotation of a genome and a genome sequence, make a pool of oligo sequences that cover all 3' UTRs.
#This can be subsetted to only consider some genes (e.g. only genes expressed in certain cells).
#Also, the parameters for oligo coverage (oligo length, coverage depth, etc.) can be tuned.

#python3

import gffutils
import os
import sys
from operator import itemgetter
from itertools import groupby
import subprocess
import pybedtools
import pickle
from Bio import SeqIO
import gzip
import argparse
import pickle
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
	#Make gff database

	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')

	polyAsites = {} #{transcriptid : [chrm, start, stop, strand]}
	genecounter = 0

	for gene in genes:
		geneshortname = gene.attributes['gene_name'][0]

		if gene.attributes['gene_type'][0] != 'protein_coding' and geneshortname != 'MALAT1':
			continue
		genecounter +=1
		if genecounter % 5000 == 0:
			print('Gene {0}...'.format(genecounter))

		for transcript in db.children(gene, featuretype = 'transcript'):
			#Only look at protein coding transcripts
			if transcript.attributes['transcript_type'][0] != 'protein_coding' and geneshortname != 'MALAT1':
				continue
			try:
				tags = transcript.attributes['tag']
				#Don't consider transcripts with these tags
				if 'cds_start_NF' in tags or 'mRNA_start_NF' in tags or 'cds_end_NF' in tags or 'mRNA_end_NF' in tags:
					continue
			except KeyError: #'tag' isn't in the attributes for this transcript
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
	command = ['liftOver', 'polyAsites.bed', '/beevol/home/taliaferro/Annotations/LiftOver/hg38ToMm10.over.chain.gz', 'liftedpolyAsites.bed', 'unliftedpolyAsites.bed']
	subprocess.call(command)
	num_lifted = sum(1 for line in open('liftedpolyAsites.bed'))
	num_unlifted = sum(1 for line in open('unliftedpolyAsites.bed')) / 2 #In unliftedpolyAsites.bed, there are 2 lines per polyAsite
	percentlifted = round((num_lifted / float(num_lifted + num_unlifted)) * 100, 2)
	print('Lifted over {0} percent of {1} polyAsites.'.format(percentlifted, num_unlifted + num_lifted))
	os.remove('polyAsites.bed')

	#Define "valid" polyAsite regions.
	#Within 200 nt of a lifted over site.  Strand does not matter.
	validwindows = {} #{chrm : [list of valid nt]}
	with open('liftedpolyAsites.bed', 'r') as infh:
		for line in infh:
			line = line.strip().split('\t')
			chrm = line[0]
			start = int(line[1])
			windowstart = start - 250 #Can change this
			windowstop = start + 250#Can change this
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
	print('Getting gene names...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)

	genes = db.features_of_type('gene')

	ens2genename = {} #{ensid : gene short name}

	for gene in genes:
		genename = str(gene.id).split('.')[0]
		shortname = gene.attributes['gene_name'][0]
		ens2genename[genename] = shortname

	print('Found names for {0} genes.'.format(len(ens2genename)))

	return ens2genename



#Iterate through a gff, gene by gene, collecting all 3' UTR coords for protein coding transcripts of that gene
def iterategff(gff, validwindows, stepsize):
	#Make gff database

	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	if os.path.isfile(db_fn) == False:
		gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	genes = db.features_of_type('gene')
	genecounter = 0
	geneswithoutpassingtranscripts = []
	transcriptcounter = 0
	longUTRcounter = 0
	tslcounter = 0
	notinvalidwindowscounter = 0

	mergedexons = {} #{genename : [chrm, strand, [[mergedexon1start, mergedexon1stop], [mergedexon2start, mergedexon2stop]], [junctionblocks]]}
	geneswewant = ['Trp53inp2', 'Trak2', 'Cplx2', 'Net1', 'Cdc42bpg', 'Gdf11', 'Malat1', 'Afap1l1', 'Ogt', 'Paxbp1', 'Fnbp4', 'Rab13', 'Akap12', 'Kif5b']

	for gene in genes:
		genename = str(gene.id).split('.')[0]
		geneshortname = gene.attributes['gene_name'][0]
		if geneshortname not in geneswewant:
			continue

		chrm = str(gene.chrom)
		if gene.attributes['gene_type'][0] != 'protein_coding' and geneshortname != 'Malat1':
			continue
		genecounter +=1
		if genecounter % 5000 == 0:
			print('Gene {0}...'.format(genecounter))

		#Does any UTR in this gene have introns
		hasjunctions = False


		#First get all coordinates that are UTR coordinate for any transcript of this gene.
		geneUTRcoords = []
		for transcript in db.children(gene, featuretype = 'transcript'):
			
			#Only look at protein coding transcripts with the exception of the one Malat1 transcript we want
			if transcript.attributes['transcript_type'][0] != 'protein_coding' and str(transcript.id) != 'ENSMUST00000172812.2':
				continue
			try:
				tsl = transcript.attributes['transcript_support_level'][0]
			except KeyError: #if the transcript doesn't have a transcript_support_level annotation
				continue
			tags = transcript.attributes['tag']
			#Don't consider transcripts with these tags
			if 'cds_start_NF' in tags or 'mRNA_start_NF' in tags or 'cds_end_NF' in tags or 'mRNA_end_NF' in tags:
				continue

			transcriptcounter +=1

			#Only consider transcripts with transcript_support_level of 1 or less
			#We really want Cplx2 and Cdc42bpg and Malat1 but their sole transcripts are tsl 2 and tsl 5 and tsl NA, respectively
			if (tsl == 'NA' or int(tsl) > 1) and (geneshortname != 'Cplx2' and geneshortname != 'Cdc42bpg' and geneshortname != 'Malat1'):
				tslcounter +=1
				#print 'TSL fail, {0}.'.format(str(transcript.id))
				continue

			#If the polyA site for this transcript is not a "valid" polyA site as defined by validwindows, forget it
			strand = transcript.strand
			if strand == '+':
				polyAsite = transcript.end
			elif strand == '-':
				polyAsite = transcript.start
			if polyAsite not in validwindows[chrm]:
				#notinvalidwindowscounter +=1
				pass
				#continue

			exoncoords = [] #[list of coords that are exonic]
			CDScoords = [] #[list of coords that are CDS exonic]

			if geneshortname != 'Malat1':
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

				#Get 100 nt immediately upstream of the beginning of the UTR.  We are going to include this as UTR sequence 
				#so that we get good coverage of the beginning of the UTR.
				#If the UTR is shorter than one oligo size (260), we need extra sequence to fit one whole oligo.
				#The amount of extra sequence we need is 110 - UTRlength
				if strand == '+':
					UTRstartindex = exoncoords.index(max(CDScoords))
					upstreamnt = exoncoords[UTRstartindex - 99 : UTRstartindex + 1]
					UTRl = (UTRend - UTRstart) + 1
					if UTRl < 260:
						extrantneeded = 260 - UTRl
						upstreamnt = exoncoords[UTRstartindex - 259 - extrantneeded : UTRstartindex + 1]
						#print UTRl, len(upstreamnt), strand 
				elif strand == '-':
					UTRstartindex = exoncoords.index(min(CDScoords))
					upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 99]
					UTRl = (UTRend - UTRstart) + 1
					if UTRl < 260:
						extrantneeded = 260 - (UTRend - UTRstart + 1)
						upstreamnt = exoncoords[UTRstartindex : UTRstartindex + 99 + extrantneeded]
						#print UTRl, len(upstreamnt), strand 

				UTRcoords = list(range(UTRstart, UTRend + 1))
				#Add on the upstream nt
				UTRcoords = upstreamnt + UTRcoords

				#Filter for those that are exonic
				UTRcoords = sorted(list(set(UTRcoords).intersection(exoncoords)))

				#Add 100 nt onto the 3' end of the UTR so that we get good coverage
				if strand == '+':
					lastnt = UTRcoords[-1]
					extrant = list(range(lastnt + 1, lastnt + 101))
					UTRcoords = UTRcoords + extrant
				elif strand == '-':
					lastnt = UTRcoords[0]
					extrant = list(range(lastnt - 100, lastnt))
					UTRcoords = extrant + UTRcoords

				#If this UTR is longer than 10 kb or shorter than 96 (160 - one step size) nt, forget it
				if len(UTRcoords) > 10000:
					print('Long UTR : {0}'.format(str(transcript.id)))
					longUTRcounter +=1
					continue

				#Now get breaks in consecutive exonic positions
				#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
				UTRexoncoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
				for k, g in groupby(enumerate(UTRcoords), lambda ix: ix[0] - ix[1]):
					exonbp = list(map(itemgetter(1), g))
					if len(exonbp) > 1:
						UTRexoncoords.append([exonbp[0], exonbp[-1]])

				#If this is Malat1, we already know what the "UTR" coords are (in this case just the transcript ends)
			elif geneshortname == 'Malat1':
				UTRexoncoords = [[transcript.start - 50, transcript.end + 50]]

			#Add this exon(s) to geneUTRcoords
			for exon in UTRexoncoords:
				geneUTRcoords.append([gene.chrom, str(exon[0]), str(exon[1])])

			#If this UTR has more than one exon, we need to make junction oligos for every exon/exon boundary
			if len(UTRexoncoords) == 1:
				pass
			elif len(UTRexoncoords) > 1:
				hasjunctions = True
				junctionblocks = [] #This is the entire region covering a junction, from 110-stepsize upstream to 110-stepsize downstream
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

					#Assuming an oligo length of 260, the most leftward junction oligo would have 260-stepsize nt on the left exon and 
					#stepsize nt on the right exon.  The most rightward junction oligo would have stepsize nt on the left exon and 260 - stepsize
					#nt on the right exon.  So what we really want is 259-stepsizeupstream---lastexonic--firstexonic--259-stepsizedownstream.
					#Then make oligos from that later.
					#If 259-stepsize away from the junction is outside the bounds of the transcript, add sequence from beyond the transcript to make up the difference.

					#Assuming an oligo length of 260, the evenly spaced bridge should be 129upstream--last exonic--first exonic--129 downstream
					#Upstream tiled should be offset by stepsize from above: 129 + stepsize upstream--lastexonic--firstexonic--129 - stepsize downstream
					#Downstream tiled should be offset by stepsize the other way: 129 - stepsize upstream--lastexonic--firstexonic--129 + stepsizedownstream
					#So what we really want is 129 + stepsize upstream--lastexonic--firstexonic--129 + stepsize downstream as a block.  Then make oligos from that later.
					#If 98 nt away from the junction is outside the bounds of the transcript, add sequence from beyond the transcript to make up the difference.

					if lastexonicindex - (259 - stepsize) < 0: #if the junction is within 259-stepsize of the transcript start
						missingntleft = ((259 - stepsize) + 1) - lastexonicindex
						junctionblockstart = exoncoords[0] - missingntleft
						junctionblock = list(range(junctionblockstart, exoncoords[0])) + exoncoords[0:lastexonicindex + 1] + exoncoords[firstexonicindex : firstexonicindex + (259 + stepsize) + 1]

					elif firstexonicindex + (259 - stepsize) > len(exoncoords) - 1: #if the junction is within 109-stepsize of the transcript end
						rightntintranscript = len(exoncoords[firstexonicindex:])
						missingntright = ((259 - stepsize) + 1) - rightntintranscript
						junctionblockend = exoncoords[-1] + missingntright
						junctionblock = exoncoords[lastexonicindex - (259 - stepsize) : lastexonicindex + 1] + exoncoords[firstexonicindex:] + list(range(exoncoords[-1] + 1, junctionblockend + 1))

					elif lastexonicindex - (259 - stepsize) >= 0 and firstexonicindex + (259 - stepsize) <= len(exoncoords) - 1: #if the entire junction block is contained within the transcript
						junctionblockstart = lastexonicindex - (259 - stepsize)
						junctionblockend = firstexonicindex + (259 - stepsize) + 1
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
					'''

					junctionblocks.append(junctionblock)


		#Write geneUTRcoords to file so we can use Bedtools
		writeUTRcoords(geneUTRcoords)

		#Merge the UTR exons
		bed = pybedtools.BedTool('temp.sorted.bed')
		mergedbed = bed.merge()

		#Get the merged "exons"
		genename = str(gene.id).split('.')[0]
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

	print('Looked through {0} protein coding genes. {1} of these ({2}%) had at least one transcript that passed all filters.'.format(genecounter, genecounter - len(geneswithoutpassingtranscripts), round(((genecounter - len(geneswithoutpassingtranscripts)) / float(genecounter)) * 100), 3))
	print('Looked through {0} transcripts. {1} had UTRs that were too long or too short. {2} were not supported in the annotation (TSL > 1). {3} did not have polyA sites within 200 nt of the orthologous human genome location.'.format(transcriptcounter, longUTRcounter, tslcounter, notinvalidwindowscounter))
	os.remove('temp.bed')
	os.remove('temp.sorted.bed')
	
	return mergedexons

def filtergenes(genesyouwant, mergedexons):
	#Taking in a file of genes you want to make oligos from, filter merged exons to only contain those genes.

	#Figure out how many genes are in this file
	genes = []
	with open(genesyouwant, 'r') as infh:
		genecounter = 0
		for line in infh:
			line = line.strip()
			genecounter +=1
			genes.append(line)

	print('Looking for {0} genes...'.format(genecounter))

	#For genes in mergedexons, if that gene isn't in expressedgenes, remove it
	for gene in list(mergedexons):
		if gene not in genes:
			del mergedexons[gene]

	print('Found {0} of {1} genes that we looked for.'.format(len(mergedexons), genecounter))
	return mergedexons


def breakoligo(coords):
	#Get breaks in consecutive oligo positions
	#http://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
	brokencoords = [] #[[exon1start, exon1stop], [exon2start, exon2stop]]
	for k, g in groupby(enumerate(coords), lambda ix: ix[0] - ix[1]):
		exonbp = list(map(itemgetter(1), g))
		if len(exonbp) > 1:
			brokencoords.append([exonbp[0], exonbp[-1]])

	return brokencoords


def cutexons(coords, stepsize):
	#For an entry in mergedexons, split it up into oligos
	strand = coords[1]
	exons = coords[2]
	junctions = coords[3]
	oligos = [] #[[oligo1start, oligo1end], [oligo2start, oligo2end]]
	junctionoligos = [] #[[j1half1start, j1half1end], [j1half2start, j1half2end], [j2half1start, j2half1end], [j2half2start, j2half2end]]

	#If there's only one exon, this is easy.
	if len(exons) >= 1:
		if strand == '+':
			for exon in exons:
				exonstart = exon[0]
				exonend = exon[1]
				currentpos = exonstart
				#This is confusing but because the gff is 1-based and open, an oligo that is x to x + 109 will actually be length 110
				while currentpos + 259 <= exonend:
					oligostart = currentpos
					oligoend = currentpos + 259
					oligos.append([[oligostart, oligoend]])
					currentpos += stepsize 
				#If the final oligo does not lie flush with the end of the exon, make one more oligo that does
				#Not really needed since we have 100nt pads
				#if currentpos != exonend:
					#oligos.append([[exonend - 259, exonend]])

		elif strand == '-':
			#Gotta start from the other end of coords so that we make the same kind of "final" oligo that lies flush with the polyA site
			#If we didn't consider strand and treated the two strands the same, the final flush oligo for the - strand would actually be flush with the stop codon
			#or 50 nt upstream of the stop codon where the exon starts
			exons = reversed(exons)
			for exon in exons:
				exonstart = exon[1] #these coords have start > end because they are on minus strand
				exonend = exon[0]
				currentpos = exonstart
				while currentpos - 259 >= exonend:
					oligoend = currentpos #these coords are back to having start < end
					oligostart = currentpos - 259
					oligos.append([[oligostart, oligoend]])
					currentpos -= stepsize 
				#If the final oligo does not lie flush with the end of the exon, make one more oligo that does
				#Not really needed since we have 100nt pads
				#if currentpos != exonend:
					#oligos.append([[exonend, exonend + 259]])

	#Wait why?  Isn't this the point of junctionblocks?  You should never have to break one of these oligos.
	#If there's more than one exon, it's more complicated
	'''
	elif len(exons) > 1:
		if strand == '+':
			exoniccoords = []
			for exon in exons:
				exoniccoords += range(exon[0], exon[1] + 1)
			currentposindex = 0
			#This is slightly different than the single-exon example above. Here we are stepping along indices, so we need the oligo to go from index x to index x + 110.
			#This will give coordinates that are 109 "places" apart
			while currentposindex + 260 <= len(exoniccoords):
				oligostart = currentposindex
				oligoend = currentposindex + 260
				oligocoords = exoniccoords[oligostart : oligoend]
				#We might need to break this oligo because it may have crossed an exon/exon boundary
				oligocoords = breakoligo(oligocoords) #if this doesn't cross a boundary, its now [[start, stop]]. if it does, its [[start1, stop1], [start2, stop2], [start3, stop3]]
				oligos.append(oligocoords) 
				currentposindex += stepsize
			#If the final oligo does not lie flush with the end of the last exon, make one more oligo that does
			if currentposindex != len(exoniccoords):
				oligocoords = exoniccoords[-260:]
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)

		elif strand == '-':
			exoniccoords = []
			for exon in exons:
				exoniccoords += range(exon[0], exon[1] + 1)
			exoniccoords.reverse()
			currentposindex = 0
			while currentposindex + 260 <= len(exoniccoords):
				oligostart = currentposindex + 260 #have to make it so start < end
				oligoend = currentposindex
				oligocoords = exoniccoords[oligoend : oligostart]
				#We might need to break this oligo because it may have crossed an exon/exon boundary
				oligocoords.reverse() #flip it back around so that startcoord < endcoord
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)
				currentposindex += stepsize
			#If the final oligo does not lie flush with the end of the last exon, make one more oligo that does
			if currentposindex != len(exoniccoords):
				oligocoords = exoniccoords[-260:]
				oligocoords = reversed(oligocoords)
				oligocoords = breakoligo(oligocoords)
				oligos.append(oligocoords)
	'''
	for junction in junctions:
		joligoindicieslist = []
		for i in list(range(len(junction)))[::stepsize]:
			joligoindicies = junction[i : i + 260]
			if len(joligoindicies) == 260:
				joligoindicieslist.append(joligoindicies)

		for oligo in joligoindicieslist:
			brokencoords = breakoligo(oligo)
			junctionoligos.append(brokencoords)

	return oligos, junctionoligos

def makeoligos(mergedexons, stepsize):
	#Go through merged exons.
	#Feed every gene to cutexons to receive back oligo coordinates
	oligocoords = {} #{gene : [chrm, strand, [oligos]]} #For oligos that do not have introns, i.e. continuous oligos
	junctionoligocoords = {} #{gene : [chrm, strand, [junction oligos]]} #For oligos that do have introns

	for gene in mergedexons:
		coords = mergedexons[gene]
		chrm = coords[0]
		strand = coords[1]
		oligos, junctionoligos = cutexons(coords, stepsize)
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
	with open('temp.gff', 'w') as outfh:
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
					elif oligolength == 259:
						if strand == '+':
							oligoend += 1
						elif strand == '-':
							oligostart -= 1
					oligocounter +=1
					oligooutlist = [chrm, 'mm10', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=regular_oneexon' + ';Parent={0}.UTR'.format(gene)]
					oligooutlists.append(oligooutlist)
				elif len(oligo) > 1: #if this oligo spans multiple exons
					oligocounter +=1
					numberofpieces = len(oligo)
					oligostart = oligo[0][0]
					oligoend = oligo[-1][1]
					oligooutlist = [chrm, 'mm10', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=regular_multiexon' + ';Parent={0}.UTR'.format(gene)]
					piececounter = 0
					oligooutlists.append(oligooutlist)
					for piece in oligo:
						piececounter +=1
						piecestart = piece[0]
						pieceend = piece[1]
						oligooutlist = [chrm, 'mm10', 'oligopiece', str(piecestart), str(pieceend), '.', strand, '.', 'ID={0}.{1}.{2}'.format(gene, oligocounter, piececounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=regularmultiexon' + ';Parent={0}.{1}'.format(gene, oligocounter)]
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
					oligooutlist = [chrm, 'mm10', 'oligo', str(oligostart), str(oligoend), '.', strand, '.', 'ID={0}.{1}'.format(gene, oligocounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=junction' + ';Parent={0}.UTR'.format(gene)]
					oligooutlists.append(oligooutlist)
					for piece in oligo:
						piececounter +=1
						piecestart = piece[0]
						pieceend = piece[1]
						oligooutlist = [chrm, 'mm10', 'junctionpiece', str(piecestart), str(pieceend), '.', strand, '.', 'ID={0}.{1}.{2}'.format(gene, oligocounter, piececounter) + ';gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';oligo_id={0}.{1}'.format(gene, oligocounter) + ';oligo_type=junction' + ';Parent={0}.{1}'.format(gene, oligocounter)]
						oligooutlists.append(oligooutlist)

			#Make gene-level line
			#This isn't actually the UTR.  It includes 50 nt or so of the coding sequence. It may not encompass the entire upstream coding region as a junction oligo piece may be further upstream than the first regular oligo start.
			geneline = [chrm, 'mm10', 'UTR', str(UTRstart), str(UTRend), '.', strand, '.', 'gene_id={0}'.format(gene) + ';gene_name={0}'.format(genenames[gene]) + ';number_of_oligos={0}'.format(oligocounter) + ';ID={0}.UTR'.format(gene)]

			outfh.write(('\t').join(geneline) + '\n')
			for oligooutlist in oligooutlists:
				outfh.write(('\t').join(oligooutlist) + '\n')


def sortoligogff(gff, outfile):
	#Because in the original gff produced by makegff the junctionoligos are added after all of the non-junction oligos have been written,
	#this puts their oligo IDs out of order with respect to genomic position.  Junction oligos account for ~5% of all oligos so it's not a big deal,
	#but still, we would like to fix this.

	#The easiest way to do this is to read in the original gff, sort the oligos for each UTR, then write a new one.

	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	utrs = db.features_of_type('UTR')

	with open(outfile, 'w') as outfh:
		for utr in utrs:
			gene = str(utr.id).split('.')[0]
			outfh.write(('\t').join([str(utr.chrom), 'mm10', 'UTR', str(utr.start), str(utr.end), '.', str(utr.strand), '.', 'gene_id={0};gene_name={1};number_of_oligos={2};ID={3}'.format(utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['number_of_oligos'][0], utr.attributes['ID'][0])]) + '\n')
			if utr.strand == '+':
				oligocounter = 0
				for oligo in db.children(utr, featuretype = 'oligo', order_by = 'start'):
					oligocounter +=1
					outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'oligo', str(oligo.start), str(oligo.end), '.', str(oligo.strand), '.', 'ID={0}.{1};gene_id={2};gene_name={3};oligo_id={4}.{5};oligo_type={6};Parent={7}'.format(utr.attributes['gene_id'][0], oligocounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], oligo.attributes['Parent'][0])]) + '\n')
					if oligo.attributes['oligo_type'][0] == 'regular_multiexon':
						piececounter = 0
						for oligopiece in db.children(oligo, featuretype = 'oligopiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'oligopiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')
					if oligo.attributes['oligo_type'][0] == 'junction':
						piececounter = 0
						for oligopiece in db.children(oligo, featuretype = 'junctionpiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'junctionpiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}.{9}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')

			elif utr.strand == '-':
				oligocounter = 0
				for oligo in db.children(utr, featuretype = 'oligo', order_by = 'end', reverse = True):
					oligocounter +=1
					outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'oligo', str(oligo.start), str(oligo.end), '.', str(oligo.strand), '.', 'ID={0}.{1};gene_id={2};gene_name={3};oligo_id={4}.{5};oligo_type={6};Parent={7}'.format(utr.attributes['gene_id'][0], oligocounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], oligo.attributes['Parent'][0])]) + '\n')
					if oligo.attributes['oligo_type'][0] == 'regular_multiexon':
						piececounter = 0
						for oligopiece in db.children(oligo, featuretype = 'oligopiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'oligopiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')
					if oligo.attributes['oligo_type'][0] == 'junction':
						piececounter = 0
						for oligopiece in db.children(oligo, featuretype = 'junctionpiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'junctionpiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}.{9}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')

	os.remove(db_fn)

def gff2fasta(gff, genomefasta):
	print('Indexing gff...')
	gff_fn = gff
	db_fn = os.path.abspath(gff_fn) + '.db'
	gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

	db = gffutils.FeatureDB(db_fn)
	print('Done indexing!')

	print('Indexing genome sequence...')
	seq_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(genomefasta, 'rt'), 'fasta'))
	print('Done indexing!')

	oligos = db.features_of_type('oligo')
	seqs = {} #{oligoname : seq}
	oligonames = [] #list of oligoIDs, this is so we can order the oligos when we write the fasta
	adapter5 = 'GCTTCGATATCCGCATGCTA'
	adapter3 = 'CTCTTGCGGTCGCACTAGTG'

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

	with open('oligos.mm10.fa', 'w') as outfh:
		for oligoname in oligonames:
			seq = seqs[oligoname]
			seq = adapter5 + seq + adapter3
			outfh.write('>' + oligoname + '\n' + seq + '\n')

	os.remove(db_fn)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff', type = str, help = 'GFF file of genome to make oligos from.')
	parser.add_argument('--othergenomegff', type = str, help = 'GFF file of genome in which we will look for orthologous polyA sites.')
	parser.add_argument('--filterdata', type = str, help = 'File containing genes to make oligos from.')
	parser.add_argument('--genomefasta', type = str, help = 'Sequence in fasta format for the genome to make oligos from.')
	parser.add_argument('--stepsize', type = int, help = 'Step size from one oligo to the next.')
	args = parser.parse_args()

	#Find orthologous polyA sites
	#validwindows = getothergenomepolyA(args.othergenomegff)
	#pickle.dump(validwindows, open('validwindows.pkl', 'wb'))
	validwindows = pickle.load(open('validwindows.pkl', 'rb'))
	#Find UTRs that pass transcript and polyA site filters
	mergedexons = iterategff(args.gff, validwindows, args.stepsize)
	pickle.dump(mergedexons, open('mergedexons.pkl', 'wb'))
	#mergedexons = pickle.load(open('mergedexons.pkl', 'rb'))
	mergedexons = filtergenes(args.filterdata, mergedexons)

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

	print('After filtering, end up with {0} nt of UTR sequence.'.format(UTRsize))

	oligocoords, junctionoligocoords = makeoligos(mergedexons, args.stepsize)

	#Print how many oligos this gets broken into
	oligocount = 0
	for gene in oligocoords:
		oligos = oligocoords[gene][2]
		oligocount += len(oligos)
	for gene in junctionoligocoords:
		oligos = junctionoligocoords[gene][2]
		oligocount += len(oligos)

	print('This is broken up into {0} oligos.'.format(oligocount))

	#Make oligo gff
	makegff(oligocoords, junctionoligocoords, genenames)
	#Sort that gff
	sortoligogff('temp.gff', 'oligos.mm10.gff')
	os.remove('temp.gff')
	#Make fasta
	gff2fasta('oligos.mm10.gff', args.genomefasta)
