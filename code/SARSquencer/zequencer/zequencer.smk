# import modules
import os
import re
import sys
from Bio import SeqIO
from Bio import Entrez
import subprocess
from shutil import copyfile
import pyfasta
import shutil
import gzip
import fileinput

# example invokation
# snakemake \
# --snakefile zequencer/zequencer.smk \
# --config \
# r1_fastq=test/50-m239-1e6_S5_L001_R1_001.fastq.gz \
# r2_fastq=test/50-m239-1e6_S5_L001_R2_001.fastq.gz \
# sample_name=50-m239-1e6 \
# ncbi_accession=M33262 \
# downsampling_fasta=test/SIVmac239-amplicons.fasta \
# reads_per_amplicon=1000 \
# reads_per_sample=10000 \
# minimum_variant_percentage=0.05 \
# minimum_read_length=100

# required config parameters
r1_fastq = config['r1_fastq'] # gzip-compressed FASTQ R1 file
r2_fastq = config['r2_fastq'] # gzip-compressed FASTQ R2 file
sample_name = config['sample_name'] # name of sample to use in output files
ncbi_accession = config['ncbi_accession'] # NCBI ncbi_accession number for reference sequence
downsampling_fasta = config['downsampling_fasta'] # (Optional) FASTA file used to map amplicon reads prior to downsampling.
reads_per_amplicon = config['reads_per_amplicon'] # (Optional) Downsample each amplicon to this number of reads
reads_per_sample = config['reads_per_sample'] # (Optional) Downsample each sample to this number of reads. Used for non-amplicon datasets primarily.
minimum_variant_percentage = config['minimum_variant_percentage'] # Minimum varaiant frequency to report. Default to 0.05.
minimum_read_length = config['minimum_read_length'] # exclude reads from mapping shorter than this value
minimum_coverage = config['minimum_coverage']

# set defaults for parameters that are passed without a value
if minimum_variant_percentage == '' : minimum_variant_percentage = '0.05'
if minimum_read_length == '' :  minimum_read_length = '100'
if reads_per_sample == '' : reads_per_sample = '100000000000000'

	
rule all:
    input:
        'out/' + sample_name + '.bam',
        'out/' + sample_name + '.vcf',
        'out/' + ncbi_accession + '.gbk',
    run:
        # remove temporary files
        shutil.rmtree('tmp/')
        shutil.rmtree('ref/')
        subprocess.call('cp out/{0}.bam ./{0}.bam'.format(sample_name),shell=True)
        subprocess.call('cp out/{0}.vcf ./{0}.vcf'.format(sample_name),shell=True)
        subprocess.call('mv out {0}'.format(sample_name),shell=True)
        subprocess.call('tar -czvf {0}.tar.gz {0}'.format(sample_name),shell=True)


## Skip the reference  getting if the files exist

rule extract_reference_files:
	message: 
		"""extract reference files"""
	input:
		ncbi_accession + '.tar.gz'
	output:
		'ref/' + ncbi_accession + '.config',
		'ref/' + ncbi_accession +'/snpEffectPredictor.bin',
		'ref/' + ncbi_accession +'/genes.gbk',   
		'ref/' + ncbi_accession + '.fa',
		'ref/' + ncbi_accession + '.fa.fai',
		'ref/' + ncbi_accession + '.fa.nix',
		'ref/' + ncbi_accession + '.gbk'
	run:
		print('untar file')
		subprocess.call('tar -xvf {0}'.format(input[0]), shell=True)
		print('rename .fasta to .fa or as needed')
		if os.path.isfile('ref/{0}.fasta'.format(ncbi_accession)):
			subprocess.call('mv ref/{0}.fasta ref/{0}.fa'.format(ncbi_accession), shell=True)
		if os.path.isfile('ref/{0}.fasta.nix'.format(ncbi_accession)):
			subprocess.call('mv ref/{0}.fasta.nix ref/{0}.fa.nix'.format(ncbi_accession), shell=True)
		if os.path.isfile('ref/{0}.fasta.fai'.format(ncbi_accession)):
			subprocess.call('mv ref/{0}.fasta.fai ref/{0}.fa.fai'.format(ncbi_accession), shell=True)
		touch(output[0])
		touch(output[1])
		touch(output[2])
		touch(output[3])
		touch(output[4])
		touch(output[5])

if not os.path.isfile('{0}.tar.gz'.format(ncbi_accession)):
	rule download_ncbi_reference:
		message: 
			"""download reference sequence from NCBI in Genbank format"""
		output:
			'ref/' + ncbi_accession + '.gbk',
			temp('ref/tmp.gbk'),
			temp('ref/tmp_cleaned.gbk'),
			'ref/' + ncbi_accession + '.fa'
		run:
			def getGenbank():
				'''retrieve Genbank file from NCBI and copy to temporary file'''

				Entrez.email = 'dhoconno@wisc.edu'  # Always tell NCBI who you are

				# Downloading...
				net_handle = Entrez.efetch(db="nucleotide", id=ncbi_accession, rettype="gb", retmode="text")
				out_handle = open(output[1], "w")
				out_handle.write(net_handle.read())
				out_handle.close()
				net_handle.close()
				print("Saved")

				# with at least some sequences (e.g., KU501215) the NCBI sequence record contains a versioning .1 or .2suffix
				# this is treated inconsistently by some tools
				# to eliminate this issue, remove versioning suffix from temporary genbank file
				# in testing, discovered that snpEff uses chromsome name from LOCUS field of Genbank file
				# other tools use ncbi_accession field for chromosome name
				# use a regular expression to change the value of the LOCUS field to match value of ncbi_accession field
				# biopython expects an exact number of spaces between fields, so we need to calculate the number of spaces to add after the replaced ncbi_accession number

				# there is almost certainly a more elegant way to do this with biopython

				with open(output[1]) as infile, open(output[2], 'w') as outfile:
					for line in infile:
						# remove version info
						line = re.sub(ncbi_accession + '\.[1-9]', ncbi_accession, line)

						# overwrite locus field with ncbi_accession number

						# number of spaces after ncbi_accession
						space_ct = 23 - len(ncbi_accession)
						spacer = ' ' * space_ct

						line = re.sub('LOCUS(\s*)\S*\s*', 'LOCUS' + r'\1' + ncbi_accession + spacer, line)
						outfile.write(line)

			def createReferenceGenbankFiles():
				'''
				Create Genbank and FASTA versions of reference genomes
				'''

				# read Genbank file and extract sequence ID
				seq_id = SeqIO.read(output[2], "genbank").id

				# create Genbank file
				SeqIO.convert(output[2], "genbank", output[0], "genbank")

				# create FASTA file
				SeqIO.convert(output[2], "genbank", output[3], "fasta")

			getGenbank()
			createReferenceGenbankFiles()

	rule index_ncbi_reference:
		message:
			"""run samtools faidx to create FASTA index"""
		input:
			'ref/' + ncbi_accession + '.fa'
		output:
			'ref/' + ncbi_accession + '.fa.fai'
		run:
			subprocess.call(['samtools', 'faidx', input[0]])

	rule create_novoindex:
		message:
			"""create Novoalign index"""
		input:
			'ref/' + ncbi_accession + '.fa'
		output:
			'ref/' + ncbi_accession + '.fa.nix'
		run:
			subprocess.call(['novoindex',output[0],input[0]])

	rule create_snpeff_database:
		message:
			"""create database for functione variant annotation with snpEff"""
		input:
			'ref/' + ncbi_accession + '.gbk',
		output:
			'ref/' + ncbi_accession + '.config',
			'ref/' + ncbi_accession +'/snpEffectPredictor.bin',
			'ref/' + ncbi_accession +'/genes.gbk',
		run:
			# get sequence id and description from Genbank file
			record = SeqIO.read(input[0], "genbank")
			seq_id = record.id
			seq_description = record.description

			# write snpEff config
			with open(output[0], "w") as myfile:
				myfile.write('data.dir = .\n')
				myfile.write('lof.ignoreProteinCodingAfter  : 0.95\n')
				myfile.write('lof.ignoreProteinCodingBefore : 0.05\n')
				myfile.write('lof.deleteProteinCodingBases : 0.50\n')
				myfile.write('codon.Standard                                                          : TTT/F, TTC/F, TTA/L, TTG/L+, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/L+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G\n')
				myfile.write('\n' + seq_id + '.genome : ' + seq_description)

			# create temporary Genbank file named genes.gbk as needed by snpeff
			copyfile(input[0], output[2])

			# create snpeff database
			cmd = ['snpEff', 'build', '-c', output[0], '-genbank', seq_id]

			subprocess.call(cmd)

			print(' '.join(cmd))

# if reads_per_amplicon and downsampling_fasta are empty then simply merge R1 and R2 FASTQ to prepare for mapping
# if reads_per_amplicon and downsampling_fasta are specified, run normalize_amplicon_coverage rule
if reads_per_amplicon == '' and downsampling_fasta == '':
    rule simple_merge_fastq:
        message:
            """run reformat.sh to create a single FASTQ file from paired-end inputs"""
        input:
            r1_fastq,
            r2_fastq
        output:
            temp('out/merged.fq.gz'),
            temp('tmp/merged.fq.gz')
        run:
            left_trim = '22'
            if not os.path.exists('tmp'): os.makedirs('tmp')
            # run bbmerge
            subprocess.call(['bbmerge.sh',
                                'in=' + input[0],
                                'in2=' + input[1],
                                'forcetrimleft=' + left_trim,
                                'out=' + output[1]])
            subprocess.call(['reformat.sh',
								'in=' + output[1],
								'out=' + output[0],
								'minlength=' + '317'])
            
            print('--Read normalizing (w/o downsampling) is 100% complete.--')
else:
    rule normalize_amplicon_coverage:
        message:
            """run normalize amplicon coverage script and output merged.fq.gz file that has been downsampled"""
        input:
            r1_fastq,
            r2_fastq,
            downsampling_fasta
        output:
            temp('out/merged.fq.gz'),
            temp('tmp/merged.fq.gz')
        run:
            # for virus sequences generated with Nick Loman ZIKV technique using multiplexed pools with many small amplicons
            # normalizes coverage for each amplicon by mapping full readset to small amplicon and retaining only subset of mapped reads
            # this is essential for normalizing coverage differences in the 30 amplicons that comprise one genome
            # original version of this tool in exp 18442
            # this version eliminates mapping of normalized reads to reference because exp 18582 workflow for mapping and variant calling is better

            # create temporary file directories if they do not already exist
            tmp_dir = 'tmp/'

            # create temporary directory and subdirectories
            if not os.path.exists('tmp'): os.makedirs('tmp')
            # if not os.path.exists(tmp_dir + '/split_reference_fasta'): os.makedirs(tmp_dir + '/split_reference_fasta')
            # if not os.path.exists(tmp_dir + '/mapped_reads'): os.makedirs(tmp_dir + '/mapped_reads')
            # if not os.path.exists(tmp_dir + '/filtered_reads'): os.makedirs(tmp_dir + '/filtered_reads')

            left_trim = '22'

            # run bbmerge
            subprocess.call(['bbmerge.sh',
                                'in=' + input[0],
                                'in2=' + input[1],
                                'forcetrimleft=' + left_trim,
                                'out=' + output[1]])
            subprocess.call(['reformat.sh',
								'in=' + output[1],
								'out=' + output[0],
								'minlength=' + '317'])
            
            print('--Read normalizing (w/o downsampling) is 100% complete.--')

rule preprocess_merged_reads:
    message:
        """remove low quality, short reads, and downsample merged.fq.gz to specified number of reads.
        This is useful when you have a whole genome dataset that has more coverage than can be supported by amount of input template"""
    input:
        'out/merged.fq.gz'
    output:
        'tmp/preprocessed.fq',
    run:
        # use bbmap reformat.sh to remove low quality sequences and prune short sequences
        # can also downsample FASTQ to appropriate number of reads
        # will also output intereleaved, decompressed FASTQ regardless of input FASTQ
        # this streamlines subsequent steps by ensuring input is uncompressed FASTQ

        print('--Preprocessing FASTQ reads with bbmap reformat.sh--')

        reformat_cmd = ['reformat.sh',
                    'in=' + input[0],
                    'out=' + output[0],
                    'qtrim=t',
                    'minlength=' + str(minimum_read_length),
                    'samplereadstarget=' + str(reads_per_sample),
                    'sampleseed=3']

        subprocess.call(reformat_cmd)

rule map_reads_to_reference:
    message:
        """Map preprocessed reads to reference with Novoalign"""
    input:
        'tmp/preprocessed.fq',
        'ref/' + ncbi_accession + '.fa.nix',
        'ref/' + ncbi_accession + '.fa.fai',
        'ref/' + ncbi_accession + '.fa'
    output:
        'tmp/mapping.filtered.sam'
    run:
    	subprocess.call('bbmap.sh in=' + input[0] + ' \
    	outm=' + output[0] + ' \
    	ref=' + input[3], shell=True)
		

rule call_variants:
    message:
        """sort and index SAM file, then use bbmap callvariants.sh to call variants"""
    input:
        'tmp/mapping.filtered.sam',
        'ref/' + ncbi_accession + '.fa',
    output:
        'out/' + sample_name + '.bam',
        'tmp/' + sample_name + '.unannotated.vcf',
    run:
        # sort SAM file
        print('--Sort SAM file and convert to BAM file--') 
        with open(output[0], "wb") as out:
            subprocess.call(['samtools',
                            'sort',
                            input[0],
                            '--reference',
                            input[1]], stdout=out)
        with open(output[1], "wb") as out: 
            variantscmd = ['callvariants.sh',
                                'in=' + output[0],
                                'ref=' + input[1],
                                'minallelefraction=' + str(minimum_variant_percentage),
                                'rarity=' + str(minimum_variant_percentage),
                                'coverage=t',
                                'calldel=t',
                                'callins=t',
                                'callsub=t',
                                'mincov=' + str(minimum_coverage),
                                'minreads=10',
                                'minvarcopies=1',
                                'vcf=' + output[1],
                                'minscore=10',
                                'overwrite=t']
       	    subprocess.call(variantscmd, stdout=out)

rule annotate_variants:
    message:
        """annotate VCF variants with snpEff"""
    input:
        'tmp/' + sample_name + '.unannotated.vcf',
        'ref/' + ncbi_accession + '.config'
    output:
        'out/' + sample_name + '.vcf'
    run:
        # run snpEff to annotate vcf file
        # only annotate variants within features (by setting ud = 0)

        with open(output[0], "wb") as out:

            snpeff_cmd = ['snpEff',
                        '-c',
                        input[1],
                        '-ud',
                        '-onlyProtein',
                        ncbi_accession,
                        input[0]]

            print(' '.join(snpeff_cmd))

            subprocess.call(snpeff_cmd, stdout=out)

        # replace generic Sample1 sample identifier with sample_name
        with fileinput.FileInput(output[0], inplace=True) as file:
            for line in file:
                print(line.replace('Sample1', sample_name), end='')


        # remove snpEff summary file and text summary invokation location
        os.remove('snpEff_summary.html')
        os.remove('snpEff_genes.txt')

rule copy_reference_genbank:
    message:
        """copies reference Genbank file to output folder so it can be loaded in Geneious"""
    input:
        'ref/' + ncbi_accession + '.gbk',
    output:
        'out/' + ncbi_accession + '.gbk',
    run:
        copyfile(input[0], output[0])