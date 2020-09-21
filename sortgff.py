import gffutils
import os
import sys

def sortoligogff(gff, outfile):
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
						for oligopiece in db.children(oligo, featuretype = 'oligopiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'junctionpiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')

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
						for oligopiece in db.children(oligo, featuretype = 'oligopiece', order_by = 'start'):
							piececounter +=1
							outfh.write(('\t').join([str(oligo.chrom), 'mm10', 'junctionpiece', str(oligopiece.start), str(oligopiece.end), '.', str(oligopiece.strand), '.', 'ID={0}.{1}.{2};gene_id={3};gene_name={4};oligo_id={5}.{6};oligo_type{7};Parent={8}'.format(utr.attributes['gene_id'][0], oligocounter, piececounter, utr.attributes['gene_id'][0], utr.attributes['gene_name'][0], utr.attributes['gene_id'][0], oligocounter, oligo.attributes['oligo_type'][0], utr.attributes['gene_id'][0], oligocounter)]) + '\n')

	os.remove(db_fn)

sortoligogff(sys.argv[1], sys.argv[2])
