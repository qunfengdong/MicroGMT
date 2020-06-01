#!/usr/bin/python

import time
import os
import sys
import subprocess
import argparse
from functions import *

def main():
	parser = argparse.ArgumentParser(description='Vcf file annotation', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-i', type=str, dest='in_dir', \
		required=True, \
		help='Input directory')
	group1.add_argument('-o', type=str, dest='out_dir', \
		required=True, \
		help='Output directory')

	group2 = parser.add_argument_group('Optional arguments')
	group2.add_argument('-r', type=str, dest='reference', \
		default='NC_045512', \
		help='Name of reference. [NC_045512]')
	group2.add_argument('-p', type=str, dest='table_prefix', \
		default='out_summary', \
		help='Prefix of summary tables for annotated vcf files. Do not include path, except for folder name(s) inside output directory!')
	group2.add_argument('-c', dest='csvS', \
		action='store_true', \
		help='Create CSV format snpEff summary files.')
	group2.add_argument('-l', type=str, dest='log', default="Annotate_vcf.log", \
		help='Name of the log file [Annotate_vcf.log]')
	group2.add_argument('-f', dest='table_form', 
		default="both", choices=['a','b','both'], \
		help='Format of summary tables [a]')
	group2.add_argument('-rg', type=str, dest='region_file', \
		default=None, \
		help='Name of the region file (Optional)')
	group2.add_argument('-na', dest='skip_anno', \
		action='store_true', \
		help='Skip vcf annotation step, just make summary tables from annotated vcfs.')
	group2.add_argument('-eff', type=str, dest='path_to_snpEff', \
		default=None, \
		help='Absolute path to snpEff.jar. Required if annotae vcf files.')

	args = parser.parse_args()

	param = {}
	param['skip_anno'] = args.skip_anno
	if args.csvS:
		param['csvS'] = "T"
	else:
		param['csvS'] = "F"
	param['in_dir'] = os.path.join(os.getcwd(),args.in_dir)
	param['out_dir'] = os.path.join(os.getcwd(),args.out_dir)
	param['path_to_snpEff'] = args.path_to_snpEff
	param['reference'] = args.reference
	param['log_nm'] = args.log
	param['out_log']=os.path.join(param['out_dir'],param['log_nm'])
	param['table_prefix'] = args.table_prefix
	param['table_form'] = args.table_form
	
	if args.region_file:
		param['reg_file']=os.path.join(os.getcwd(),args.region_file)
	else:
		param['reg_file']=None

	with open(param['out_log'],'w') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] Annotate_vcf started.\n')

	f.close()

	log_print(param['out_log'],'==================== MicroGMT ====================')
	log_print(param['out_log'],'                  Annotate_vcf')
	log_print(param['out_log'],'             Version 1.2  (May 2020)')
	log_print(param['out_log'],'   Bug report: Yue Xing <yue.july.xing@gmail.com>')
	log_print(param['out_log'],'======================================================')

	if not param['skip_anno']:
		if not param['path_to_snpEff']:
			scn_print("Error: Path to snpEff not provided!")
			log_print(param['out_log'],"Error: Path to snpEff not provided!")
			exit(1)

	if not param['skip_anno']:
		log_print(param['out_log'],"Start annotating vcf files...")
		vcf_annotate(param['out_log'],param['in_dir'],param['out_dir'], \
			param['path_to_snpEff'],param['reference'],param['csvS'])
		log_print(param['out_log'],"Annotation of vcf files successful!")
		if param['table_form']=='a':
			make_out_table_form1(param['out_log'],param['out_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
		elif param['table_form']=='b':
			make_out_table_form2(param['out_log'],param['out_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
		else:
			make_out_table_form1(param['out_log'],param['out_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
			make_out_table_form2(param['out_log'],param['out_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
	else:
		log_print(param['out_log'],"Skipped annotation of vcf files.")
		if param['table_form']=='a':
			make_out_table_form1(param['out_log'],param['in_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
		elif param['table_form']=='b':
			make_out_table_form2(param['out_log'],param['in_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
		else:
			make_out_table_form1(param['out_log'],param['in_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
			make_out_table_form2(param['out_log'],param['in_dir'], \
				param['out_dir'],param['table_prefix'],param['reg_file'])
	log_print(param['out_log'],"Finished making summary tables!")

if __name__ == '__main__':
	main()



