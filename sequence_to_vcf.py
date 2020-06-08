#!/usr/bin/python

import time
import os
import sys
import subprocess
import argparse
from functions import *

def main():
	parser = argparse.ArgumentParser(description='Sequence file(s) to vcf file(s)', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-r', type=str, dest='ref_genome', required=True, \
		help='Fasta formatted reference genome file')
	group1.add_argument('-i', dest='input_format', \
		choices=['assembly','contig','fastq'], \
		required=True, help='Type of input file')
	group1.add_argument('-o', type=str, dest='out_dir', required=True, \
		help='Output directory')

	group2 = parser.add_argument_group('Additional arguments for inputs')
	group2.add_argument('-fs', type=str, dest='fasta_seqs', default=None, \
		help='Fasta format assembly or contig file.')
	group2.add_argument('-fq1', type=str, dest='fastq1', default=None, \
		help='Fastq file 1. For paired-end fastq data.')
	group2.add_argument('-fq2', type=str, dest='fastq2', default=None, \
		help='Fastq file 2. For paired-end fastq data.')
	group2.add_argument('-fq', type=str, dest='fastq', default=None, \
		help='Fastq file. For single-end fastq data.')

	group3 = parser.add_argument_group('Optional arguments')
	group3.add_argument('-l', type=str, dest='log', default="Sequence_to_vcf.log", \
		help='Name of the log file [Sequence_to_vcf.log]')
	group3.add_argument('-n', type=str, dest='name', default="test", \
		help="Name of the input sample. Does not work with 'assembly' option. [test]")
	group3.add_argument('-gatk', type=str, dest='path_to_gatk', default=None, \
		help="Absolute path to GenomeAnalysisTK.jar. Only required for 'fastq' option.")
	group3.add_argument('-picard', type=str, dest='path_to_picard', default=None, \
		help="Absolute path to picard.jar. Only required for 'fastq' option.")
	group3.add_argument('-kb', dest='keep_bam', action='store_true', \
		help='Keep BAM files.')
	group3.add_argument('-ki', dest='keep_idx', action='store_true', \
		help="Keep index files. Only works with 'fastq' option.")
	group3.add_argument('-p', type=float, dest='prior', \
		default=0, help="Prior for bcftools variant caller (expected substitution rate). 0 means the prior is disabled. Only works for 'assembly' or 'contig' option. [0]'.")
	group3.add_argument('-m', type=int, dest='mbq', 
		default=10, help="Minimum base quality for variant caller. Only works with 'fastq' option. [10]")
	group3.add_argument('-a', dest='asm', choices=['asm5','asm10','asm20'], 
		default='asm5', help="Sequence divergence: asm5/asm10/asm20 for ~0.1/1/5 percentages. Only works with 'assembly' option. [asm5]")
	group3.add_argument('-t', dest='thread', type=int, 
		default=10, help="Number of threads. [10]")

	args = parser.parse_args()

	param = {}
	param['out_dir'] = os.path.join(os.getcwd(),args.out_dir)
	param['log_nm'] = args.log
	param['out_log']=os.path.join(param['out_dir'],param['log_nm'])

	with open(param['out_log'],'w') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] Program started.\n')

	f.close()

	log_print(param['out_log'],'==================== MicroGMT ====================')
	log_print(param['out_log'],'                 Sequence_to_vcf')
	log_print(param['out_log'],'             Version 1.3  (June 2020)')
	log_print(param['out_log'],'   Bug report: Yue Xing <yue.july.xing@gmail.com>')
	log_print(param['out_log'],'======================================================')

	param['ref_genome']=os.path.join(os.getcwd(),args.ref_genome)
	param['input_format'] = args.input_format
	param['fasta_seqs'] = args.fasta_seqs
	param['fastq1'] = args.fastq1
	param['fastq2'] = args.fastq2
	param['fastq'] = args.fastq
	param['path_to_gatk'] = args.path_to_gatk
	param['path_to_picard'] = args.path_to_picard

	if param['input_format']=='assembly' \
	or param['input_format']=='contig':
		if not param['fasta_seqs']:
			scn_print("Error: Assembly or contig file not provided!")
			log_print(param['out_log'],"Error: Assembly or contig file not provided!")
			exit(1)
		else:
			param['fasta_seqs']=os.path.join(os.getcwd(),param['fasta_seqs'])
	else: 
		if not param['path_to_gatk']:
			scn_print("Error: Path to GATK not provided!")
			log_print(param['out_log'],"Error: Path to GATK not provided!")
			exit(1)
		if not param['path_to_picard']:
			scn_print("Error: Path to picard not provided!")
			log_print(param['out_log'],"Error: Path to picard not provided!")
			exit(1)
		if param['fastq'] and (not param['fastq1']) \
		and (not param['fastq2']):
			param['pairing']='single'
			param['fastq']=os.path.join(os.getcwd(),param['fastq'])
			param['fastq1']="nofile"
			param['fastq1']="nofile"
			log_print(param['out_log'],"Fastq file is single-end.")
		elif param['fastq1'] and param['fastq2'] \
		and not param['fastq']:
			param['pairing']='paired'
			param['fastq1']=os.path.join(os.getcwd(),param['fastq1'])
			param['fastq2']=os.path.join(os.getcwd(),param['fastq2'])
			param['fastq']="nofile"
			log_print(param['out_log'],"Fastq files are paired-end.")
		else:
			scn_print("Error: Fastq file(s) not provided correctly!")
			log_print(param['out_log'],"Error: Fastq file(s) not provided correctly!")
			exit(1)

	param['name'] = args.name
	param['thread'] = str(args.thread)
	
	if args.keep_bam:
		param['keep_bam'] = "T"
		log_print(param['out_log'],"BAM files will be kept.")
	else:
		param['keep_bam'] = "F"

	if args.keep_idx:
		param['keep_idx'] = "T"
		log_print(param['out_log'],"Index files will be kept for 'fastq' option.")
	else:
		param['keep_idx'] = "F"

	param['prior'] = str(args.prior)
	param['mbq'] = str(args.mbq)
	param['BAQ'] = 'not_used'
	param['asm'] = args.asm

	if param['input_format']=='assembly':
		# fasta to vcf
		# The ID in fasta header will be the name
		log_print(param['out_log'],"Start processing fasta assembly inputs.")
		fasta_to_vcf(param['out_log'],param['out_dir'], \
			param['fasta_seqs'],param['ref_genome'],param['asm'], \
			param['prior'],param['keep_bam'],param['thread'])
		log_print(param['out_log'],"Processing fasta assembly inputs successful!")
	elif param['input_format']=='contig':
		# contigs to vcf
		# The ID in fasta header will be provided by user
		log_print(param['out_log'],"Start processing fasta contig inputs.")
		contig_to_vcf(param['out_log'],param['out_dir'], \
			param['fasta_seqs'],param['ref_genome'], \
			param['name'],param['prior'],param['keep_bam'],param['thread'])
		#log_print(param['out_log'],"Finished processing fasta contig inputs.")
	elif param['input_format']=='fastq':
		# fastq to vcf
		# The ID in fasta header will be provided by user
		log_print(param['out_log'],"Start processing fastq inputs.")
		fastaq_to_vcf(param['out_log'],param['out_dir'],param['ref_genome'], \
			param['fastq1'],param['fastq2'],param['name'],param['prior'], \
			param['mbq'],param['BAQ'],param['path_to_picard'], \
			param['path_to_gatk'],param['keep_bam'], \
			param['keep_idx'],param['fastq'],param['pairing'],param['thread'])
		log_print(param['out_log'],"Processing fastq inputs successful!")

if __name__ == '__main__':
	main()
