#Given the results of oligo enrichments, as calculated by DESeq2, define windows
#of enriched oligos that encompass active RNA elements.  Then, once you have
#windows, define minimal sequences within those windows that represent the
#minimal active element. These would be sequences that are in common to every
#oligo within the window.

#Windows are defined by looking for an enriched oligo and extending a window if either
#(A) the next oligo is sig or
#(B) there are x sig oligos in a row after a gap of x oligos

#After defining these primary windows, two windows are merged together if 
#the distance between them is less than  (width of the narrow window + 2)
#After merging, any window less than 10 oligos long is removed

import os
import pandas as pd
import sys
import gffutils
from itertools import groupby
from operator import itemgetter

def simpleroligos2genome(gff):
    #Just make a table that relates oligo position in the UTR (i.e. 1, 2, 3, 4, etc.)
    #to the start genomecoord of the oligo

    oligocoords = {} #{utr : {oligoposition : genomecoordofstart}}

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
        if utrid not in oligocoords:
            oligocoords[utrid] = {}
        for oligo in db.children(utr, featuretype = 'oligo'):
            oligoid = int(str(oligo.id).split('.')[1])
            if oligo.strand == '+':
                genomecoord = int(oligo.start)
            elif oligo.strand == '-':
                genomecoord = int(oligo.stop)
            oligocoords[utrid][oligoid] = genomecoord
            if utrid == 'Afap1l1' and oligoid == 1:
                print('yo')

    #turn into df
    df = pd.concat({k: pd.DataFrame.from_dict(v, 'index') for k, v in oligocoords.items()}, axis = 0)
    #get indexes out as column names
    df.reset_index(level = 1, inplace = True)
    df.reset_index(level = 0, inplace = True)
    df.columns = ['genename', 'oligopos', 'genomecoord']

    return df


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
		genomecoords = utrcoords[utr]
		oligocoords = list(range(len(genomecoords) + 1))
		oligocoords2genomecoords = dict(zip(oligocoords, genomecoords))
		utr_oligocoords[utr] = oligocoords2genomecoords

	
	#turn utr_oligocoords into a df and return the df
	df = pd.concat({k: pd.DataFrame.from_dict(v, 'index') for k, v in utr_oligocoords.items()}, axis = 0)
	#get indexes out as column names
	df.reset_index(level = 1, inplace = True)
	df.reset_index(level = 0, inplace = True)
	df.columns = ['genename', 'utrcoord', 'genomecoord']

	return df


def definewindows(deseqtable, windowlengthminimum):
    enrichedwindows = {} #{gene : [(windowstart, windowstop)]}
    #The input for this function is a modified DESeq2 output table
    #It contains the following columns:
    #ensid, genename, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, sig, neuritelog2FC, expected
    #sig is either 'soma', 'neurite', or 'no'
    
    #Read file
    df = pd.read_csv(deseqtable, sep = '\t', header = 0, index_col = None)
    #ensid is of the form ENSMUSG00000.0
    #the number after the dot is the position of the oligo along the UTR
    df[['ensid', 'oligopos']] = df['ensid'].str.split('.', n = 1, expand = True)
    #Reorder columns
    df = df[['ensid', 'oligopos', 'genename', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'sig', 'neuritelog2FC', 'expected']]
    #Get rid of columns we don't need
    df = df[['genename', 'oligopos', 'sig']]
    #Convert oligopos to integer
    df = df.astype({'oligopos' : int})
    #Add in the start nucleotide position (0-based) of each oligo relative to the beginning of the first oligo for this gene
    #Oligos are 260 nt long, 4 nt step size
    oligopos = df.oligopos.tolist()
    utrcoord = [(x - 1) * 4 for x in oligopos]
    df = df.assign(utrcoord = utrcoord)
    #Arrange by genename then oligopos
    df = df.sort_values(by = ['genename', 'oligopos'], axis = 0)


    #OK, now for each gene, define windows of significance
    #Walk along UTR. oligopos = x. If oligo at x sig == 'neurite' oligo, open window
    #If x + 1 oligo is sig, extend window and repeat.
    #If not, go until you find a sig. Distance between current pos and next sig is y
    #If -all- of the oligos between x + y and x + 2y are sig, then extend window and go back to the beginning.
    #If not, window ends.
    #If window length passes threshold, keep it.

    genes = df.genename.tolist()
    genes = list(set(genes))
    genedfs = [] #list of individual genedfs
    for gene in genes:
        genedf = df.loc[df['genename'] == gene]
        alloligopos = genedf.oligopos.tolist()
        sigoligopos = genedf.loc[genedf['sig'] == 'neurite'].oligopos.tolist()

        inwindow = False
        oligoposinwindow = [] #oligo positions in -any- window
        currentidx = 0
        for idx, oligopos in enumerate(alloligopos):
            #print(gene, idx, oligopos, inwindow, oligoposinwindow)
            if idx < currentidx:
                continue
            if inwindow == False:
                if oligopos not in sigoligopos:
                    continue
                elif oligopos in sigoligopos:
                    inwindow = True
                    oligoposinwindow.append(oligopos)
            elif inwindow == True:
                if oligopos in sigoligopos: #if we are in a window and the next oligo is also sig, just extend the window
                    oligoposinwindow.append(oligopos)
                    continue
                elif oligopos not in sigoligopos:
                    findasigoligoaftergap = False #are we able to find a sig oligo within 10 oligos of gap start
                    for i in range(1, 11):
                        doesthisgappass = False #does this gap fit the requirements of x sig oligos after a gap of x
                        try:
                            oligoinquestion = alloligopos[idx + i]
                            if oligoinquestion in sigoligopos:# if we found a gap that is < 10
                                findasigoligoaftergap = True
                                gapendidx = alloligopos.index(oligoinquestion)
                                #All oligos from gapidx up to (but not including) gapidx + i must be sig in order for window to continue
                                postgapsigoligos = [] #booleans for the sig-ness of the required post-gap oligos
                                doesthisgappass = False #is this gap valid? Does it have the required number of consecutive sig oligos after the gap?
                                for j in range(0, i):
                                    try:
                                        x = alloligopos[gapendidx + j]
                                        if x in sigoligopos:
                                            postgapsigoligos.append(True)
                                        else:
                                            postgapsigoligos.append(False)
                                    except IndexError:
                                        postgapsigoligos.append(False)
                                if all(postgapsigoligos): #are all items in postgapsigoligos True
                                    doesthisgappass = True
                                    break
                                else:
                                    doesthisgappass = False

                        except IndexError:
                            inwindow = False
                            doesthisgappass = False
                            break

                        if doesthisgappass == True:
                            currentidx = idx + i + i
                            #Append all oligos from idx up to (but not including) currentidx to oligoposinwindow
                            for k in range(idx, currentidx):
                                oligoposinwindow.append(alloligopos[k])
                            inwindow = True
                            break
                    
                    #If we didn't find a sig oligo within 10 positions of the gapstart
                    if findasigoligoaftergap == False:
                        inwindow = False

        #Turn this list of positions in a window into a list of tuples of windowstart/windowstop
        #https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        windows = [] #tuples of (windowstart, windowstop)
        for k, g in groupby(enumerate(oligoposinwindow), lambda x:x[0] - x[1]):
            group = list(map(itemgetter(1), g))
            windows.append((group[0], group[-1]))

        #Merge windowA and windowB if they are less than min(windowAlength, windowBlength) + 2 apart
        mergedwindows = [] #list of positions in merged windows
        for windowA in windows:
            windowAstart = windowA[0]
            windowAend = windowA[1]
            windowAlength = (windowAend - windowAstart) + 1
            wasmerged = False
            for windowB in windows:
                windowBstart = windowB[0]
                windowBend = windowB[1]
                windowBlength = (windowBend - windowBstart) + 1
                allowabledistance = min([windowAlength, windowBlength]) + 2
                if abs(windowBstart - windowAend) <= allowabledistance or abs(windowAstart - windowBend) <= allowabledistance:
                    mergedwindow = (min(windowAstart, windowBstart), max(windowAend, windowBend))
                    for i in range(mergedwindow[0], mergedwindow[1] + 1):
                        mergedwindows.append(i)
                    wasmerged = True

            #if this window was not merged with anything, just put the original window into mergedwindows
            if not wasmerged:
                for i in range(windowAstart, windowAend + 1):
                    mergedwindows.append(i)

        #remove duplicates and sort mergedwindows
        mergedwindows = sorted(list(set(mergedwindows)))
        
        #Break consective windows into chunks as above
        windows = [] #tuples of (windowstart, windowstop)
        for k, g in groupby(enumerate(mergedwindows), lambda x:x[0] - x[1]):
            group = list(map(itemgetter(1), g))
            windows.append((group[0], group[-1]))
        
        #Filter for windows that are at least windowlengthminimum long
        windows = [window for window in windows if ((window[1] - window[0]) + 1) >= windowlengthminimum]
        
        #Add data about windows to dictionary
        positionsinawindow = []
        enrichedwindows[gene] = []
        for window in windows:
            enrichedwindows[gene].append(window)
            for position in range(window[0], window[1] + 1):
                positionsinawindow.append(position)

        #Now add a new column to the df that tells whether or not that position is in a window
        inwindow = []
        for i in alloligopos:
            if i in positionsinawindow:
                inwindow.append('yes')
            else:
                inwindow.append('no')

        #Add column to df
        genedf = genedf.assign(inwindow = inwindow)
        genedfs.append(genedf)

    #Concatentate all genedfs back together
    allgenedfs = pd.concat(genedfs, axis = 0)
    
    return enrichedwindows, allgenedfs

 
def getminimalseqs(enrichedwindows, oligogff):
    #Given a set of oligo windows provided by definewindows(), get the minimal sequence
    #within those windows associated with localization. If windowmode == 'intersect', these will be defined as the sequence
    #in common to all oligos within the window. If windowmode == 'union', these will be defined as the sequence
    #in any oligo within the window

    #enrichedwindows has the following structure:
    #{genename : [(first oligo in window1, last oligo in window1), (first oligo in window2, last oligo in window2)]}
    #In this dictionary, oligos are defined by their position within the UTR, starting with 1
    #This designations correspond to the numbers after the ENS id and '.' in the oligo gff

    minimalseqs = [] #[[chrm, 'mm10', 'minimalseq', start, stop, '.', strand, attributes]]
    windowoligosgff = []

    print('Indexing gff...')
    gff_fn = oligogff
    db_fn = os.path.abspath(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, merge_strategy = 'merge', verbose = True)

    db = gffutils.FeatureDB(db_fn)
    print('Done indexing!')


    for gene in enrichedwindows:
        windowcounter = 0
        for window in enrichedwindows[gene]:
            windowcounter +=1
            coverednt = [] #[[all nt covered by oligo1], [all nt covered by oligo2], ...]
            firstoligo = window[0]
            lastoligo = window[1]
            utrs = db.features_of_type('UTR')
            for utr in utrs:
                genename = utr.attributes['gene_name'][0]
                if genename == gene:
                    chrm = str(utr.chrom)
                    strand = str(utr.strand)
                    ensid = str(utr.id).split('.')[0]
                    for oligo in db.children(utr, featuretype = 'oligo', level = 1):
                        oligonumber = int(str(oligo.id).split('.')[1])
                        #if this is one of the oligos we are after
                        if oligonumber >= firstoligo and oligonumber <= lastoligo:
                            #is this a oneexon oligo?
                            if oligo.attributes['oligo_type'][0] == 'regular_oneexon':
                                coverednt.append(list(range(oligo.start, oligo.end + 1)))
                                windowoligosgff.append([str(oligo.chrom), 'sigWindowOligos', 'oligo', str(oligo.start), str(oligo.stop), '.', oligo.strand, '.', 'ID={0};gene_name={1};oligo_type={2}'.format(str(oligo.id), str(oligo.attributes['gene_name'][0]), 'regular_oneexon')])
                            elif oligo.attributes['oligo_type'][0] == 'junction':
                                thisoligosnt = []
                                windowoligosgff.append([str(oligo.chrom), 'sigWindowOligos', 'oligo', str(oligo.start), str(oligo.stop), '.', oligo.strand, '.', 'ID={0};gene_name={1};oligo_type={2}'.format(str(oligo.id), str(oligo.attributes['gene_name'][0]), 'junction')])
                                for junctionpiece in db.children(oligo, featuretype = 'junctionpiece'):
                                    windowoligosgff.append([str(oligo.chrom), 'sigWindowOligos', 'junctionpiece', str(junctionpiece.start), str(junctionpiece.stop), '.', oligo.strand, '.', 'ID={0};gene_name={1};oligo_type={2};Parent={3}'.format(str(junctionpiece.id), str(junctionpiece.attributes['gene_name'][0]), 'junction', str(junctionpiece.attributes['Parent'][0]))])
                                    junctionpiecent = list(range(junctionpiece.start, junctionpiece.end + 1))
                                    thisoligosnt += junctionpiecent
                                coverednt.append(thisoligosnt)

            #now get the nt in the window
            if gene == 'Cplx2':
                windowmode = 'union'
            else:
                windowmode = 'intersect'
            if windowmode == 'intersect':
                commonnt = set(coverednt[0])
                for x in coverednt[1:]:
                    commonnt.intersection_update(x)
                commonnt = sorted(list(commonnt))
            elif windowmode == 'union':
                setlist = [set(x) for x in coverednt]
                unionnt = set().union(*setlist)
                commonnt = sorted(list(unionnt))

            #Does this minimal element span an intron?
            #Turn this list of positions in a window into a list of tuples of windowstart/windowstop
            #https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
            #length of minimalelementwindows is the number of exons contained in this minimal element
            minimalelementwindows = [] #tuples of (windowstart, windowstop)
            for k, g in groupby(enumerate(commonnt), lambda x:x[0] - x[1]):
                group = list(map(itemgetter(1), g))
                minimalelementwindows.append((group[0], group[-1]))

            #if this whole minimal element is contained within one exon
            if len(minimalelementwindows) == 1:
                minimalseq = [chrm, 'mm10', 'minimalseq', str(minimalelementwindows[0][0]), str(minimalelementwindows[-1][1]), '.', strand, '.', 'ID={0}.minimalseq{1};gene_name={2};element_type={3}'.format(ensid, windowcounter, gene, 'singleexon')]
                minimalseqs.append(minimalseq)
            #if it spans multiple exons
            elif len(minimalelementwindows) > 1:
                minimalelementexoncounter = 0
                wholewindow = [chrm, 'mm10', 'minimalseq', str(minimalelementwindows[0][0]), str(minimalelementwindows[-1][1]), '.', strand, '.', 'ID={0}.minimalseq{1};gene_name={2};element_type={3}'.format(ensid, windowcounter, gene, 'multiexon')]
                minimalseqs.append(wholewindow)
                for x in minimalelementwindows:
                    minimalelementexoncounter +=1
                    minimalseq = [chrm, 'mm10', 'minimalseqexon', str(x[0]), str(x[1]), '.', strand, '.', 'ID={0}.minimalseq{1}.exon{2};gene_name={3};element_type={4};Parent={5}'.format(ensid, windowcounter, minimalelementexoncounter, gene, 'multiexon', ensid + '.minimalseq' + windowcounter)]
                    minimalseqs.append(minimalseq)

    return minimalseqs, windowoligosgff

                

#Relate oligo coords to genomecoords
oligocoorddf = simpleroligos2genome(sys.argv[1])
print(oligocoorddf.head())


#Define enriched windows of oligos
enrichedwindows, df = definewindows(sys.argv[2], 10)
print(df.head())

#Join oligocoord and window dfs
df = pd.merge(oligocoorddf, df, how = 'inner', on = ['genename', 'oligopos'])
print(enrichedwindows)
df.to_csv(path_or_buf = os.path.abspath(sys.argv[2]) + '.sigoligowindows', sep = '\t', na_rep = 'NA', index = False, header = True)

#Define seqs of interest based on enriched windows
minimalseqs, windowoligosgff = getminimalseqs(enrichedwindows, sys.argv[1])
with open('minimalseqs.gff', 'w') as outfh:
    for minimalseq in minimalseqs:
        print(int(minimalseq[4]) - int(minimalseq[3]))
        outfh.write(('\t').join(minimalseq) + '\n')

with open('windowoligos.gff', 'w') as outfh:
    for windowoligo in windowoligosgff:
        outfh.write(('\t').join(windowoligo) + '\n')

#TODO