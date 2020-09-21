#We've decided that we would like to make pools that tile really densely.  But a potential problem is that 
#reads that we get from this screen may be really difficult to map uniquely to oligos from these pools.
#To explore this, make a pseudo oligo file with a given step size, then make a pseudo fastq of sequences 
#made from that oligo file.  You will know how many sequences came from each file.  Then map with bowtie2
#and see how well it does.

import random
import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pysam


def makeoligos(testseq, noligos, stepsize):
	noligos = int(noligos)
	stepsize = int(stepsize)

	#Make an oligo fasta file from testseq (testseq.fa is human chromosome 1 with the Ns removed)
	with open(testseq, 'r') as infh, open('simoligos.fa', 'w') as outfh:
		for line in infh:
			line = line.strip()
			if '>' in line:
				continue
			else:
				giantseq = line

		#Sizes of the fake UTRs that our fake oligos are coming from
		#Not a great simulation of actual UTRs
		utrsizes = range(500, 5000)
		
		chosenoligos = 0
		UTRnumber = 0
		while chosenoligos < noligos:
			UTRnumber +=1
			#Pick a utr size
			utrsize = random.choice(utrsizes)
			#Pick a start site
			possiblestartsites = range(len(giantseq) - utrsize)
			startsite = random.choice(possiblestartsites)
			utr = giantseq[startsite : startsite + utrsize]

			#now make the oligos from this utr
			oligoseqs = []
			oligostart = 0
			while oligostart + 110 < len(utr):
				oligoseq = utr[oligostart : oligostart + 110] #110-mer
				oligoseqs.append(oligoseq)
				oligostart += stepsize

			#Dump oligos to file
			adapter5 = 'GCTTCGATATCCGCATGCTA'
			adapter3 = 'CTCTTGCGGTCGCACTAGTG'
			#Not adding adapters
			for i, oligo in enumerate(oligoseqs):
				outfh.write('>' + 'UTR{0}.oligo{1}'.format(UTRnumber, i) + '\n' + oligo + '\n')

			chosenoligos += len(oligoseqs)


def decideabundances(fasta):
	#Given a fasta of "oligos" decide on abundances for each of them

	#Get number of oligos
	numberofoligos = 0
	with open(fasta, 'r') as infh:
		for record in SeqIO.parse(fasta, 'fasta'):
			numberofoligos +=1

	randomnumbers = []
	for i in range(numberofoligos):
		randomnumbers.append(random.choice(range(10000)))

	#Turn random numbers into fractions (number / sum of numbers)

	abundances = [randomnumber / sum(randomnumbers) for randomnumber in randomnumbers]

	#Assign these abundances to oligos
	oligoabundances = {} #{oligo : abundance}
	with open(fasta, 'r') as infh:
		oligonumber = 0
		for record in SeqIO.parse(fasta, 'fasta'):
			oligoabundances[str(record.id)] = abundances[oligonumber]
			oligonumber +=1

	return oligoabundances


def makequalityscore(seq):
	#given a sequence, make a random quality score string (all scores >= 30)
	#possibleqscores = ['?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
	possibleqscores = range(30, 41)
	qscores = random.choices(possibleqscores, k = len(seq))

	return qscores

def makedeletions(seq, delrate):
	#given a sequence, make 1nt deletions at a rate of delrate per nt
	positionstodelete = []
	for i in range(len(seq)):
		shouldidelete = np.random.choice(['yes', 'no'], p = [delrate, 1-delrate])
		if shouldidelete == 'yes':
			positionstodelete.append(i)

	if not positionstodelete:
		return seq
	else:
		s = bytearray(seq, 'utf8')
		for ind, position in enumerate(positionstodelete):
			del s[position - ind] #position-ind is because when something is deleted, everyone's indices change

		return(s.decode('utf8'))

def makemutations(seq, mutrate):
	#given a sequence, make mutations at a rate of mutrate per nt
	mutdict = {'A' : ['C', 'T', 'G'], 'T' : ['A', 'C', 'G'], 'C' : ['A', 'G', 'T'], 'G' : ['T', 'C', 'A']}

	positionstomutate = []
	for i in range(len(seq)):
		shouldimutate = np.random.choice(['yes', 'no'], p = [mutrate, 1-mutrate])
		if shouldimutate == 'yes':
			positionstomutate.append(i)

	if not positionstomutate:
		return seq
	else:
		s = bytearray(seq, 'utf8')
		for position in positionstomutate:
			currentnt = chr(s[position])
			mutantnt = random.choice(mutdict[currentnt])
			s[position] = ord(mutantnt)

		return(s.decode('utf8'))


def makefastq(oligosequences, oligoabundances, numberofreads, mutrate, delrate):
	#Make a simulated fastq of the sequences in fasta based on the abundances in oligoabundances
	#error rate is the per nt error rate
	#del rate in the per nt deletion rate (a deletion of 1 nt)

	#We are going to assume that deletions are SYNTHESIS errors and therefore model them on the oligo prior to make reads from it
	#We are going to assume that mutations are SEQUENCING errors and therefore model them on the read

	actualabundances = {} #{oligo id : number of reads made from it}

	with open(oligosequences, 'r') as infh:
		forrecords = []
		revrecords = []
		readcounter = 0
		for record in SeqIO.parse(infh, 'fasta'):
			oligoseq = str(record.seq)
			#Make deletions
			oligoseq = makedeletions(oligoseq, delrate)

			numberofreadstomake = round(oligoabundances[str(record.id)] * numberofreads)
			actualabundances[str(record.id)] = numberofreadstomake

			for i in range(numberofreadstomake):
				readcounter += 1
				if readcounter % 100000 == 0:
					print(readcounter)
				
				#Make forward read record
				forreadseq = oligoseq[:97]
				forreadseq = makemutations(forreadseq, mutrate)
				qualityscores = makequalityscore(forreadseq)
				forrecord = SeqRecord(Seq(forreadseq, generic_dna), id = 'Read{0}.R1_{1}.read{2}'.format(readcounter, str(record.id), i), description = '')
				forrecord.letter_annotations['phred_quality'] = qualityscores
				forrecords.append(forrecord)

				#Make reverse read record
				revreadseq = oligoseq[-91:]
				revreadseq = makemutations(revreadseq, mutrate)
				qualityscores = makequalityscore(revreadseq)
				revrecord = SeqRecord(Seq(revreadseq, generic_dna).reverse_complement(), id = 'Read{0}.R2_{1}.read{2}'.format(readcounter, str(record.id), i), description = '')
				revrecord.letter_annotations['phred_quality'] = qualityscores
				revrecords.append(revrecord)

	#Shuffle the lists of records, BUT KEEP THE ORDER THE SAME BETWEEN THEM BECAUSE THEY ARE PAIRED
	allreads = list(zip(forrecords, revrecords))
	random.shuffle(allreads)
	forrecords, revrecords = zip(*allreads)

	SeqIO.write(forrecords, 'simfastq.R1.fastq', 'fastq')
	SeqIO.write(revrecords, 'simfastq.R2.fastq', 'fastq')

	with open('actualabundances.txt', 'w') as outfh:
		for oligo in actualabundances:
			outfh.write(('\t').join([oligo, str(actualabundances[oligo])]) + '\n')

def analyzesam(sam):
	#After mapping the reads (with bowtie2), see how often the read is mapped to the correct oligo.
	#This can be done easily by just looking at the sam.  The oligo that the read was made from is indicated in the read name.
	#So we can ask if that oligo matches the reference (oligo) that it was mapped to.

	with pysam.AlignmentFile(sam, 'r') as infh:
		readcounter = 0
		correctlymappedreads = 0
		abundances = {} #{oligoname : [actual reads, mapped reads]}
		mapqs = {} #{'correctlymapped' : [list of mapqs from correctly mapped reads], 'incorrectlymapped' : [list of mapqs from incorrectly mapped reads]}
		mapqs['correctlymapped'] = []
		mapqs['incorrectlymapped'] = []

		for read in infh.fetch(until_eof = True):
			if read.is_read1: #only forward reads
				readcounter +=1
				if readcounter % 1000000 == 0:
					print('Read {0}...'.format(readcounter))

				correctoligo = read.query_name.split('_')[1].split('.')
				correctoligo = ('.').join(correctoligo[0:2])
				mappedoligo = read.reference_name
				mapq = int(read.mapping_quality)
				if mappedoligo == correctoligo:
					correctlymappedreads +=1
					mapqs['correctlymapped'].append(mapq)
				elif mappedoligo != correctoligo:
					mapqs['incorrectlymapped'].append(mapq)
				if correctoligo not in abundances:
					abundances[correctoligo] = [0, 0]
				if mappedoligo not in abundances:
					abundances[mappedoligo] = [0, 0]
				abundances[correctoligo][0] +=1
				abundances[mappedoligo][1] +=1

	print('{0} of {1} reads were mapped correctly ({2}%).'.format(correctlymappedreads, readcounter, round((correctlymappedreads / readcounter) * 100, 2)))

	with open('samresults.txt', 'w') as outfh:
		outfh.write(('\t').join(['oligo', 'actualcounts', 'mappedcounts']) + '\n')
		for oligo in abundances:
			outfh.write(('\t').join([oligo, str(abundances[oligo][0]), str(abundances[oligo][1])]) + '\n')
	with open('mapqs.txt', 'w') as outfh:
		outfh.write(('\t').join(['readstatus', 'mapq']) + '\n')
		for mapq in mapqs['correctlymapped']:
			outfh.write(('\t').join(['correctlymapped', str(mapq)]) + '\n')
		for mapq in mapqs['incorrectlymapped']:
			outfh.write(('\t').join(['incorrectlymapped', str(mapq)]) + '\n')			





#Make a fastq (10M reads)
#set an error rate (0.002)
#set an indel rate (0.0001)
#Might as well make it a 'post-trimming' fastq


#makeoligos(sys.argv[1], sys.argv[2], sys.argv[3])
#oligoabundances = decideabundances('simoligos.fa')
#makefastq('simoligos.fa', oligoabundances, 10000000, 0.002, 0.0001)

analyzesam(sys.argv[1])