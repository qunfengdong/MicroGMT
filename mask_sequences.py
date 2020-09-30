#!/usr/bin/python

import time
import os
import sys
import subprocess
import argparse
from functions import *

def main():
	parser = argparse.ArgumentParser(description='Mask sequences', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-i', type=str, dest='in_table', \
		required=True, \
		help='Input summary table (Only ".all.form2.txt" tables)')
	group1.add_argument('-m', type=str, dest='mask', \
		required=True, \
		help='Input mask file')
	group1.add_argument('-d', type=str, dest='out_dir', \
		required=True, \
		help='Output directory')

	group2 = parser.add_argument_group('Optional arguments')
	group2.add_argument('-p', type=str, dest='out_pref', \
		default='masked', \
		help='Prefix of the output annotated summary tables. Do not include path, except for folder name(s) inside output directory!')
	group2.add_argument('-l', type=str, dest='log', default="Mask_sequences.log", \
		help='Name of the log file [Mask_sequences.log]')

	args = parser.parse_args()

	param = {}
	param['out_dir'] = os.path.join(os.getcwd(),args.out_dir)
	param['out_pref'] = args.out_pref
	param['in_table'] = os.path.join(os.getcwd(),args.in_table)
	param['in_mask'] = os.path.join(os.getcwd(),args.mask)
	param['out_log']=os.path.join(param['out_dir'],args.log)

	with open(param['out_log'],'w') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] Mask sequences started.\n')

	f.close()

	log_print(param['out_log'],'==================== MicroGMT ====================')
	log_print(param['out_log'],'                   Mask_sequences')
	log_print(param['out_log'],'             Version 1.3  (June 2020)')
	log_print(param['out_log'],'   Bug report: Yue Xing <yue.july.xing@gmail.com>')
	log_print(param['out_log'],'======================================================')

	mask_seqs(param['in_table'],param['in_mask'],param['out_dir'],param['out_pref'])
	log_print(param['out_log'],"Successfully completed!")

if __name__ == '__main__':
	main()
