#!/usr/bin/python

import time
import os
import sys
import subprocess
import argparse
from functions import *

def main():
	parser = argparse.ArgumentParser(description='Remove strains from table', \
		formatter_class=argparse.RawTextHelpFormatter)
	group1 = parser.add_argument_group('Mandatory inputs')
	group1.add_argument('-i', type=str, dest='in_table', \
		required=True, \
		help='Input summary table (Only ".all.form1.txt" or ".all.form2.txt" table is required!)')
	group1.add_argument('-r', type=str, dest='remove_list', \
		required=True, \
		help='List of IDs to be removed from table')
	group1.add_argument('-d', type=str, dest='out_dir', \
		required=True, \
		help='Output directory')
	group1.add_argument('-f', dest='form', \
		required=True, choices=['a','b'], \
		help='Format of summary tables')

	group2 = parser.add_argument_group('Optional arguments')
	group2.add_argument('-p', type=str, dest='out_pref', \
		default='removed', \
		help='Prefix of the output summary tables.Do not include path, except for folder name(s) inside output directory! ')
	group2.add_argument('-l', type=str, dest='log', default="Remove_from_summary_tables.log", \
		help='Name of the log file [Remove_from_summary_tables.log]')

	args = parser.parse_args()

	param = {}
	param['in_table'] = os.path.join(os.getcwd(),args.in_table)
	param['remove_list'] = os.path.join(os.getcwd(),args.remove_list)
	param['out_dir'] = os.path.join(os.getcwd(),args.out_dir)
	param['form'] = args.form
	param['out_pref'] = args.out_pref
	param['log_nm'] = args.log
	param['out_log']=os.path.join(param['out_dir'],param['log_nm'])

	with open(param['out_log'],'w') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] Removal of strains from summary tables started.\n')

	f.close()

	log_print(param['out_log'],'==================== MicroGMT ====================')
	log_print(param['out_log'],'            Remove_from_summary_tables')
	log_print(param['out_log'],'             Version 1.2  (May 2020)')
	log_print(param['out_log'],'   Bug report: Yue Xing <yue.july.xing@gmail.com>')
	log_print(param['out_log'],'======================================================')

	if param['form']=="a":
		log_print(param['out_log'],"Summary table form 1.")
		remove_from_table_form1(param['out_dir'],param['out_pref'], \
			param['in_table'],param['remove_list'])

	if param['form']=="b":
		log_print(param['out_log'],"Summary table form 2.")
		remove_from_table_form2(param['out_dir'],param['out_pref'], \
			param['in_table'],param['remove_list'])

	log_print(param['out_log'],"Finished removing strains from summary tables!")

if __name__ == '__main__':
	main()







