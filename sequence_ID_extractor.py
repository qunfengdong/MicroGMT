#!/usr/bin/python

import time
import os
import sys
import subprocess
import argparse
from functions import *

def main():
	parser = argparse.ArgumentParser(description='sequence ID extractor', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-i', type=str, dest='in_table', \
		required=True, \
		help='Input summary table (format 2 table)')
	group1.add_argument('-o', type=str, dest='out_table', \
		required=True, \
		help='Processed output table')
	group1.add_argument('-id', type=str, dest='seq_id', \
		required=True, \
		help='Strain/sequence ID to extract')
	group1.add_argument('-f', dest='out_form', \
		required=True, choices=['l','s'], \
		help='The output table format (l: long format, s: short format)')
	group1.add_argument('-a', dest='include_annot', \
		required=True, choices=['y','n'], \
		help='The summary table include custom annotation or not? (y: yes, n: no)')


	group2 = parser.add_argument_group('Optional arguments')
	group2.add_argument('-l', type=str, dest='log', default="ID_extraction.log", \
		help='Directory and name of the log file [ID_extraction.log]')

	args = parser.parse_args()

	param = {}
	param['in_table'] = os.path.join(os.getcwd(),args.in_table)
	param['out_table'] = os.path.join(os.getcwd(),args.out_table)
	param['fid'] = args.seq_id
	param['out_log'] = args.log
	param['out_form'] = args.out_form
	param['include_annot'] = args.include_annot

	with open(param['out_log'],'w') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] Strain/sequence ID extraction started.\n')

	f.close()

	log_print(param['out_log'],'==================== MicroGMT ====================')
	log_print(param['out_log'],'              sequence_ID_extractor')
	log_print(param['out_log'],'             Version 1.3.2  (Sep 2020)')
	log_print(param['out_log'],'   Bug report: Yue Xing <yue.july.xing@gmail.com>')
	log_print(param['out_log'],'======================================================')

	if param['include_annot']=="n":
		mutation_summary(param['out_log'],param['in_table'],param['fid'],param['out_table'],param['out_form'])
	else:
		mutation_summary_wanno(param['out_log'],param['in_table'],param['fid'],param['out_table'],param['out_form'])

	log_print(param['out_log'],"Successfully completed!")

if __name__ == '__main__':
	main()
