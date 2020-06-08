#!/usr/bin/python
# MicroGMT Version 1.3  (June 2020)

import time
import os
import sys
import subprocess
import argparse

def scn_print(message):
	print('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] ' + message + "\n")
	sys.stdout.flush()

def log_print(log,message):
	with open(log,'a') as f:
		f.write('[MicroGMT ' + \
			time.asctime( time.localtime(time.time()) ) + \
			'] ' + message + "\n")
	f.close()

def fasta_to_vcf(out_log,out_dir,in_seqs_file,in_ref, \
	asm,prior,keep_bam,in_th):
	in_id_list_file="id.list"
	in_id_list=os.path.join(out_dir,in_id_list_file)
	in_seqs=in_seqs_file.split("/")[-1]
	dirn = os.getcwd()
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	subprocess.call(['chmod', 'u=rwx', 'Prepare_fasta.sh'])
	subprocess.call(['./Prepare_fasta.sh',in_ref,in_seqs_file,\
		in_seqs,out_dir,asm,prior,out_log,in_th,keep_bam])
	os.chdir(dirn)

def contig_to_vcf(out_log,out_dir,in_seqs_file,in_ref, \
	in_fid,prior,keep_bam,in_th):
	dirn = os.getcwd()
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	subprocess.call(['chmod', 'u=rwx', 'Prepare_contigs.sh'])
	subprocess.call(['./Prepare_contigs.sh',in_ref,in_seqs_file,\
		in_fid,out_dir,prior,out_log,in_th])
	os.chdir(dirn)
	if keep_bam=="F":
		subprocess.call(['rm', '-f',os.path.join(out_dir,in_fid+'.bam')])

def fastaq_to_vcf(out_log,out_dir,in_ref, \
	in_fq1,in_fq2,in_fid,prior,mbq,BAQ, \
	path_to_picard,path_to_gatk,keep_bam,keep_idx, \
	in_fastq,in_pairing,in_th):
	picard_path = os.path.join(path_to_picard,"picard.jar")
	ref_name=in_ref.split("/")[-1]
	ref_name=ref_name.split(".")[0:-1]
	ref_name='.'.join(ref_name)
	ref_genome=os.path.join(out_dir,ref_name+'.fa')
	ref_genome_dict = os.path.join(out_dir, ref_name+'.dict')
	subprocess.call(['cp',in_ref,ref_genome])
	log_print(out_log,'Create index files for the reference...')
	if not os.path.exists(ref_genome_dict):
		subprocess.call(['java', '-jar', picard_path, 'CreateSequenceDictionary', \
			'REFERENCE=' + ref_genome, 'OUTPUT=' + ref_genome_dict])
	if not os.path.exists(os.path.join(out_dir, ref_name+'.fa.sa')):
		subprocess.call(['bwa', 'index', ref_genome])
	if not os.path.exists(os.path.join(out_dir, ref_name+'.fa.fai')):
		subprocess.call(['samtools','faidx',ref_genome])
	if not os.path.exists(ref_genome_dict):
		log_print(out_log,'Error: Fail to create index (.dict) for the reference genome.')
		exit(1)
	if not os.path.exists(os.path.join(out_dir, ref_name+'.fa.fai')):
		log_print(out_log,'Error: Fail to create index (.fai) for the reference genome.')
		exit(1)
	if not os.path.exists(os.path.join(out_dir, ref_name+'.fa.sa')):
		log_print(out_log,'Error: Fail to create bwa indexes for the reference genome.')
		exit(1)
	dirn = os.getcwd()
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	subprocess.call(['chmod', 'u=rwx', 'Prepare_fastq.sh'])
	subprocess.call(['./Prepare_fastq.sh',ref_genome,\
		in_fq1,in_fq2,in_fid,prior,mbq,BAQ,\
		path_to_gatk,out_dir,out_log, \
		in_fastq,in_pairing,in_th,os.path.dirname(os.path.realpath(__file__))])
	os.chdir(dirn)
	if keep_bam=="F":
		subprocess.call(['rm', '-f',os.path.join(out_dir,in_fid+'.bam')])
		subprocess.call(['rm', '-f',os.path.join(out_dir,in_fid+'.bam.bai')])
	if keep_idx=="F":
		subprocess.call(['rm', '-f',ref_genome_dict])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.fai')])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.sa')])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.pac')])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.bwt')])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.ann')])
		subprocess.call(['rm', '-f',os.path.join(out_dir, ref_name+'.fa.amb')])

def vcf_annotate(out_log,in_dir,out_dir,snpeff_dir,in_db_ref,in_csvS):
	in_tdir=os.path.dirname(os.path.realpath(__file__))
	dirn = os.getcwd()
	os.chdir(in_tdir)
	subprocess.call(['chmod', 'u=rwx', 'Annotate_vcf.sh'])
	subprocess.call(['./Annotate_vcf.sh', \
		in_dir,out_dir,snpeff_dir,in_db_ref,in_tdir,out_log,in_csvS])
	os.chdir(dirn)

def write_tabs_form1(fun_tab,fun_tab_out_table,fun_hd):
	names=list(fun_tab.keys())
	with open(fun_tab_out_table,'w') as f:
		f.write("Chr\tPos\t")
		f.write("\t".join(fun_hd))
		f.write("\n")
		for name in names:
			f.write("\t".join(name) + "\t")
			f.write("\t".join(fun_tab[name]))
			f.write("\n")
	f.close()

def write_tabs_form2(fun_tab,fun_tab_out_table,fun_hd):
	names=list(fun_tab.keys())
	with open(fun_tab_out_table,'w') as f:
		f.write("Chr\tPos\tGid\tGname\t")
		f.write("\t".join(fun_hd))
		f.write("\n")
		for name in names:
			f.write("\t".join(name) + "\t")
			f.write("\t".join(fun_tab[name]))
			f.write("\n")
	f.close()

def make_out_table_form1(out_log,in_dir,out_dir,out_pref,reg_file):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form1.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form1.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form1.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form1.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form1.txt")
	tabge_out_table=os.path.join(out_dir,out_pref+".gene.form1.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form1.txt")
	tabgmut2_out_table=os.path.join(out_dir,out_pref+".gene_name_mut.form1.txt")
	tabcg={}
	tabpg={}
	tabv={}
	tabef={}
	tabge={}
	tab={}
	tabgmut={}
	tabgmut2={}
	fl=os.listdir(in_dir)
	vfls=[j for i,j in enumerate(fl) if 'anno.vcf' in j]
	idx=list(range(len(vfls)))
	for s in idx:
		vfl=vfls[s]
		#s_nm=vfl.split('.anno.vcf')[0]
		vfl=os.path.join(in_dir,vfl)
		with open(vfl) as tmp:
			for line in tmp:
				if not line.startswith("#"):
					chr,pos,id,ref,alt,qual,ft,info,fmt,sname=line.rstrip().split("\t")
					info=info.split(";")[1:]
					info=';'.join(info)
					anno = dict(i.split("=") for i in info.split(";"))["ANN"].split(",")[0].split("|")
					vari=anno[1]
					ef=anno[2]
					gname=anno[3]
					gid=anno[4]
					gby=anno[7]
					cg=anno[9].replace('c.','')
					pg=anno[10].replace('p.','')
					wpos=(chr,pos)
					if wpos not in tabge.keys():
						tabge[wpos]= ["R"] * len(idx)
					if wpos not in tabcg.keys():
						tabcg[wpos]= ["R"] * len(idx)
					if wpos not in tabpg.keys():
						tabpg[wpos]= ["R"] * len(idx)
					if wpos not in tabv.keys():
						tabv[wpos]= ["R"] * len(idx)
					if wpos not in tabef.keys():
						tabef[wpos]= ["R"] * len(idx)
					if wpos not in tab.keys():
						tab[wpos]= ["R"] * len(idx)
					if wpos not in tabgmut.keys():
						tabgmut[wpos]= ["R"] * len(idx)
					if wpos not in tabgmut2.keys():
						tabgmut2[wpos]= ["R"] * len(idx)
					tabge[wpos][s]=gid+"|"+gname
					tabgmut[wpos][s]=ref+"|"+alt
					tabgmut2[wpos][s]=gid+"|"+gname+"|"+ref+">"+alt
					tabcg[wpos][s]=gid+"|"+gname+"|"+cg
					tabpg[wpos][s]=gid+"|"+gname+"|"+pg
					tabv[wpos][s]=gid+"|"+gname+"|"+vari
					tabef[wpos][s]=gid+"|"+gname+"|"+ef
					tab[wpos][s]=gid+"|"+gname+"|"+ref+"|"+alt \
					+"|"+gby+"|"+vari+"|"+ef+"|"+cg+"|"+pg
		tmp.close()
	hd = [sub.replace(".anno.vcf","") for sub in vfls] 
	if reg_file:
		reg={}
		with open(reg_file) as tmp2:
			for line in tmp2:
				fid,region=line.rstrip().split("\t")
				reg[fid]=region
		hd2 = []
		for ele in hd:
			if ele in reg.keys():
				ele2=ele+"|"+reg[ele]
			else:
				ele2=ele+"|UnKn"
			hd2.append(ele2)
	else:
		hd2=hd
	write_tabs_form1(tabcg,tabcg_out_table,hd2)
	write_tabs_form1(tabpg,tabpg_out_table,hd2)
	write_tabs_form1(tabv,tabv_out_table,hd2)
	write_tabs_form1(tabef,tabef_out_table,hd2)
	write_tabs_form1(tab,tab_out_table,hd2)
	write_tabs_form1(tabge,tabge_out_table,hd2)
	write_tabs_form1(tabgmut,tabgmut_out_table,hd2)
	write_tabs_form1(tabgmut2,tabgmut2_out_table,hd2)

def make_out_table_form2(out_log,in_dir,out_dir,out_pref,reg_file):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form2.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form2.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form2.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form2.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form2.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form2.txt")
	tabcg={}
	tabpg={}
	tabv={}
	tabef={}
	tab={}
	tabgmut={}
	fl=os.listdir(in_dir)
	vfls=[j for i,j in enumerate(fl) if 'anno.vcf' in j]
	idx=list(range(len(vfls)))
	for s in idx:
		vfl=vfls[s]
		#s_nm=vfl.split('.anno.vcf')[0]
		vfl=os.path.join(in_dir,vfl)
		with open(vfl) as tmp:
			for line in tmp:
				if not line.startswith("#"):
					chr,pos,id,ref,alt,qual,ft,info,fmt,sname=line.rstrip().split("\t")
					info=info.split(";")[1:]
					info=';'.join(info)
					anno = dict(i.split("=") for i in info.split(";"))["ANN"].split(",")[0].split("|")
					vari=anno[1]
					ef=anno[2]
					gname=anno[3]
					gid=anno[4]
					gby=anno[7]
					cg=anno[9].replace('c.','')
					pg=anno[10].replace('p.','')
					wpos=(chr,pos,gid,gname)
					if wpos not in tabcg.keys():
						tabcg[wpos]= ["R"] * len(idx)
					if wpos not in tabpg.keys():
						tabpg[wpos]= ["R"] * len(idx)
					if wpos not in tabv.keys():
						tabv[wpos]= ["R"] * len(idx)
					if wpos not in tabef.keys():
						tabef[wpos]= ["R"] * len(idx)
					if wpos not in tab.keys():
						tab[wpos]= ["R"] * len(idx)
					if wpos not in tabgmut.keys():
						tabgmut[wpos]= ["R"] * len(idx)
					tabgmut[wpos][s]=ref+"|"+alt
					tabcg[wpos][s]=cg
					tabpg[wpos][s]=pg
					tabv[wpos][s]=vari
					tabef[wpos][s]=ef
					tab[wpos][s]=ref+"|"+alt \
					+"|"+gby+"|"+vari+"|"+ef+"|"+cg+"|"+pg
		tmp.close()
	hd = [sub.replace(".anno.vcf","") for sub in vfls] 
	if reg_file:
		reg={}
		with open(reg_file) as tmp2:
			for line in tmp2:
				fid,region=line.rstrip().split("\t")
				reg[fid]=region
		hd2 = []
		for ele in hd:
			if ele in reg.keys():
				ele2=ele+"|"+reg[ele]
			else:
				ele2=ele+"|UnKn"
			hd2.append(ele2)
	else:
		hd2=hd
	write_tabs_form2(tabcg,tabcg_out_table,hd2)
	write_tabs_form2(tabpg,tabpg_out_table,hd2)
	write_tabs_form2(tabv,tabv_out_table,hd2)
	write_tabs_form2(tabef,tabef_out_table,hd2)
	write_tabs_form2(tab,tab_out_table,hd2)
	write_tabs_form2(tabgmut,tabgmut_out_table,hd2)

def combine_dicts(dic1,dic2,l1,l2):
	for k in dic2.keys():
		if k in dic1.keys():
			dic1[k]=dic1[k]+dic2[k]
		else:
			dic1[k]=(["R"]*l1)+dic2[k]
	for x in dic1.keys():
		if x not in dic2.keys():
			dic1[x]=dic1[x]+(["R"]*l2)
	return dic1

def read_in_table_form1(fl):
	tabcg={}
	tabpg={}
	tabv={}
	tabef={}
	tabge={}
	tab={}
	tabgmut={}
	tabgmut2={}
	with open(fl) as f:
		hd = f.readline().rstrip()
		hd=hd.split("\t")[2:]
		ln=len(hd)
		for line in f:
			line=line.rstrip().split("\t")
			wpos=(line[0],line[1])
			if wpos not in tabge.keys():
				tabge[wpos]= ["R"] * ln
			if wpos not in tabcg.keys():
				tabcg[wpos]= ["R"] * ln
			if wpos not in tabpg.keys():
				tabpg[wpos]= ["R"] * ln
			if wpos not in tabv.keys():
				tabv[wpos]= ["R"] * ln
			if wpos not in tabef.keys():
				tabef[wpos]= ["R"] * ln
			if wpos not in tab.keys():
				tab[wpos]= ["R"] * ln
			if wpos not in tabgmut.keys():
				tabgmut[wpos]= ["R"] * ln
			if wpos not in tabgmut2.keys():
				tabgmut2[wpos]= ["R"] * ln
			info=line[2:]
			#print(info)
			for i in range(len(info)):
				if info[i] != "R":
					#print(info[i])
					gid,gname,ref,alt,gby,vari,ef,cg,pg=info[i].split("|")
					tabge[wpos][i]=gid+"|"+gname
					tabgmut[wpos][i]=ref+"|"+alt
					tabgmut2[wpos][i]=gid+"|"+gname+"|"+ref+">"+alt
					tabcg[wpos][i]=gid+"|"+gname+"|"+cg
					tabpg[wpos][i]=gid+"|"+gname+"|"+pg
					tabv[wpos][i]=gid+"|"+gname+"|"+vari
					tabef[wpos][i]=gid+"|"+gname+"|"+ef
					tab[wpos][i]=gid+"|"+gname+"|"+ref+"|"+alt \
					+"|"+gby+"|"+vari+"|"+ef+"|"+cg+"|"+pg
	f.close()
	return ln,hd,tabge,tabgmut,tabgmut2,tabcg,tabpg,tabv,tabef,tab

def read_in_table_form2(fl):
	tabcg={}
	tabpg={}
	tabv={}
	tabef={}
	tab={}
	tabgmut={}
	with open(fl) as f:
		hd = f.readline().rstrip()
		hd=hd.split("\t")[4:]
		ln=len(hd)
		for line in f:
			line=line.rstrip().split("\t")
			wpos=(line[0],line[1],line[2],line[3])
			if wpos not in tabcg.keys():
				tabcg[wpos]= ["R"] * ln
			if wpos not in tabpg.keys():
				tabpg[wpos]= ["R"] * ln
			if wpos not in tabv.keys():
				tabv[wpos]= ["R"] * ln
			if wpos not in tabef.keys():
				tabef[wpos]= ["R"] * ln
			if wpos not in tab.keys():
				tab[wpos]= ["R"] * ln
			if wpos not in tabgmut.keys():
				tabgmut[wpos]= ["R"] * ln
			info=line[4:]
			#print(info)
			for i in range(len(info)):
				if info[i] != "R":
					ref,alt,gby,vari,ef,cg,pg=info[i].split("|")
					tabgmut[wpos][i]=ref+"|"+alt
					tabcg[wpos][i]=cg
					tabpg[wpos][i]=pg
					tabv[wpos][i]=vari
					tabef[wpos][i]=ef
					tab[wpos][i]=ref+"|"+alt \
					+"|"+gby+"|"+vari+"|"+ef+"|"+cg+"|"+pg
	f.close()
	return ln,hd,tabgmut,tabcg,tabpg,tabv,tabef,tab

def add_to_table_form1(out_log,in_table1,in_table2,out_pref,out_dir):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form1.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form1.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form1.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form1.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form1.txt")
	tabge_out_table=os.path.join(out_dir,out_pref+".gene.form1.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form1.txt")
	tabgmut2_out_table=os.path.join(out_dir,out_pref+".gene_name_mut.form1.txt")
	in_fl=os.path.join(os.getcwd(),in_table1)
	in_fl2=os.path.join(os.getcwd(),in_table2)
	ln1,hd1,tabge1,tabgmut1,tabgmut21,tabcg1,tabpg1, \
	tabv1,tabef1,tab1 = read_in_table_form1(in_fl)
	ln2,hd2,tabge2,tabgmut2,tabgmut22,tabcg2,tabpg2, \
	tabv2,tabef2,tab2 = read_in_table_form1(in_fl2)
	tabge=combine_dicts(tabge1,tabge2,ln1,ln2)
	tabgmut=combine_dicts(tabgmut1,tabgmut2,ln1,ln2)
	tabgmut2=combine_dicts(tabgmut21,tabgmut22,ln1,ln2)
	tabcg=combine_dicts(tabcg1,tabcg2,ln1,ln2)
	tabpg=combine_dicts(tabpg1,tabpg2,ln1,ln2)
	tabv=combine_dicts(tabv1,tabv2,ln1,ln2)
	tabef=combine_dicts(tabef1,tabef2,ln1,ln2)
	tab=combine_dicts(tab1,tab2,ln1,ln2)
	hd=hd1+hd2
	write_tabs_form1(tabge,tabge_out_table,hd)
	write_tabs_form1(tabgmut,tabgmut_out_table,hd)
	write_tabs_form1(tabgmut2,tabgmut2_out_table,hd)
	write_tabs_form1(tabcg,tabcg_out_table,hd)
	write_tabs_form1(tabpg,tabpg_out_table,hd)
	write_tabs_form1(tabv,tabv_out_table,hd)
	write_tabs_form1(tabef,tabef_out_table,hd)
	write_tabs_form1(tab,tab_out_table,hd)

def add_to_table_form2(out_log,in_table1,in_table2,out_pref,out_dir):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form2.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form2.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form2.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form2.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form2.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form2.txt")
	in_fl=os.path.join(os.getcwd(),in_table1)
	in_fl2=os.path.join(os.getcwd(),in_table2)
	ln1,hd1,tabgmut1,tabcg1,tabpg1, \
	tabv1,tabef1,tab1 = read_in_table_form2(in_fl)
	ln2,hd2,tabgmut2,tabcg2,tabpg2, \
	tabv2,tabef2,tab2 = read_in_table_form2(in_fl2)
	tabgmut=combine_dicts(tabgmut1,tabgmut2,ln1,ln2)
	tabcg=combine_dicts(tabcg1,tabcg2,ln1,ln2)
	tabpg=combine_dicts(tabpg1,tabpg2,ln1,ln2)
	tabv=combine_dicts(tabv1,tabv2,ln1,ln2)
	tabef=combine_dicts(tabef1,tabef2,ln1,ln2)
	tab=combine_dicts(tab1,tab2,ln1,ln2)
	hd=hd1+hd2
	write_tabs_form2(tabgmut,tabgmut_out_table,hd)
	write_tabs_form2(tabcg,tabcg_out_table,hd)
	write_tabs_form2(tabpg,tabpg_out_table,hd)
	write_tabs_form2(tabv,tabv_out_table,hd)
	write_tabs_form2(tabef,tabef_out_table,hd)
	write_tabs_form2(tab,tab_out_table,hd)

def del_in_dict(tt,fun_idx):
	tt_new={}
	for k in tt.keys():
		tt_new[k] = [y for x,y in enumerate(tt[k]) if x not in fun_idx]
	return tt_new

def remove_from_table_form1(out_dir,out_pref, in_fl,in_rl):
	ntabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form1.txt")
	ntabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form1.txt")
	ntabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form1.txt")
	ntabef_out_table=os.path.join(out_dir,out_pref+".effect.form1.txt")
	ntab_out_table=os.path.join(out_dir,out_pref+".all.form1.txt")
	ntabge_out_table=os.path.join(out_dir,out_pref+".gene.form1.txt")
	ntabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form1.txt")
	ntabgmut2_out_table=os.path.join(out_dir,out_pref+".gene_name_mut.form1.txt")
	ln,hd,tabge,tabgmut,tabgmut2,tabcg,tabpg,tabv,tabef,tab=read_in_table_form1(in_fl)
	rml=[]
	with open(in_rl,'r') as f:
		for line in f:
			line = line.rstrip().rstrip()
			rml.append(line)
	idx = [i for i,j in enumerate(hd) if j.split("|")[0] in rml]
	hd2 = [y for x,y in enumerate(hd) if x not in idx]
	ntabge=del_in_dict(tabge,idx)
	ntabgmut=del_in_dict(tabgmut,idx)
	ntabgmut2=del_in_dict(tabgmut2,idx)
	ntabcg=del_in_dict(tabcg,idx)
	ntabpg=del_in_dict(tabpg,idx)
	ntabv=del_in_dict(tabv,idx)
	ntabef=del_in_dict(tabef,idx)
	ntab=del_in_dict(tab,idx)
	write_tabs_form1(ntabge,ntabge_out_table,hd2)
	write_tabs_form1(ntabgmut,ntabgmut_out_table,hd2)
	write_tabs_form1(ntabgmut2,ntabgmut2_out_table,hd2)
	write_tabs_form1(ntabcg,ntabcg_out_table,hd2)
	write_tabs_form1(ntabpg,ntabpg_out_table,hd2)
	write_tabs_form1(ntabv,ntabv_out_table,hd2)
	write_tabs_form1(ntabef,ntabef_out_table,hd2)
	write_tabs_form1(ntab,ntab_out_table,hd2)

def remove_from_table_form2(out_dir,out_pref, in_fl,in_rl):
	ntabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form2.txt")
	ntabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form2.txt")
	ntabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form2.txt")
	ntabef_out_table=os.path.join(out_dir,out_pref+".effect.form2.txt")
	ntab_out_table=os.path.join(out_dir,out_pref+".all.form2.txt")
	ntabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form2.txt")
	ln,hd,tabgmut,tabcg,tabpg,tabv,tabef,tab=read_in_table_form2(in_fl)
	rml=[]
	with open(in_rl,'r') as f:
		for line in f:
			line = line.rstrip().rstrip()
			rml.append(line)
	idx = [i for i,j in enumerate(hd) if j.split("|")[0] in rml]
	hd2 = [y for x,y in enumerate(hd) if x not in idx]
	ntabgmut=del_in_dict(tabgmut,idx)
	ntabcg=del_in_dict(tabcg,idx)
	ntabpg=del_in_dict(tabpg,idx)
	ntabv=del_in_dict(tabv,idx)
	ntabef=del_in_dict(tabef,idx)
	ntab=del_in_dict(tab,idx)
	write_tabs_form2(ntabgmut,ntabgmut_out_table,hd2)
	write_tabs_form2(ntabcg,ntabcg_out_table,hd2)
	write_tabs_form2(ntabpg,ntabpg_out_table,hd2)
	write_tabs_form2(ntabv,ntabv_out_table,hd2)
	write_tabs_form2(ntabef,ntabef_out_table,hd2)
	write_tabs_form2(ntab,ntab_out_table,hd2)

def change_table_format(fl,out):
	with open(out,'w') as o:
		o.write("chr\tpos\tgene_id\tgene_name\tmutation\tID\tregion\n")
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			for line in f:
				line = line.rstrip().split("\t")
				for i in range(4,len(line)):
					if line[i]!="R" and line[i]!="":
						if len(hd[i].split("|"))>1:
							new=[line[0],line[1],line[2],line[3],line[i],hd[i].split("|")[0],hd[i].split("|")[1]]
						else:
							new=[line[0],line[1],line[2],line[3],line[i],hd[i].split("|")[0],"Unknown"]
						o.write("\t".join(new)+"\n")
	f.close()
	o.close()

def find_uniq_mutations(fl,out):
	with open(out,'w') as o:
		o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tmutation\n")
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			for line in f:
				line = line.rstrip().split("\t")
				new=line[4:len(line)]
				a=[i for i,j in enumerate(new) if j!="R"]
				#b=[j for i,j in enumerate(new) if j!="R"] 
				if len(a)==1:
					s=a[0]+4
					ss=line[0:4]
					ss.append(line[s])
					if len(hd[s].split("|"))>1:
						o.write(hd[s].split("|")[0]+"\t"+hd[s].split("|")[1]+"\t"+"\t".join(ss)+"\n")
					else:
						o.write(hd[s].split("|")[0]+"\t"+"Unknown"+"\t"+"\t".join(ss)+"\n")
	f.close()
	o.close()

def change_table_format_wanno(fl,out):
	with open(out,'w') as o:
		o.write("chr\tpos\tgene_id\tgene_name\tcustom_annot\tmutation\tID\tregion\n")
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			for line in f:
				line = line.rstrip().split("\t")
				for i in range(5,len(line)):
					if line[i]!="R" and line[i]!="":
						if len(hd[i].split("|"))>1:
							new=[line[0],line[1],line[2],line[3],line[4],line[i],hd[i].split("|")[0],hd[i].split("|")[1]]
						else:
							new=[line[0],line[1],line[2],line[3],line[4],line[i],hd[i].split("|")[0],"Unknown"]
						o.write("\t".join(new)+"\n")
	f.close()
	o.close()

def find_uniq_mutations_wanno(fl,out):
	with open(out,'w') as o:
		o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tcustom_annot\tmutation\n")
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			for line in f:
				line = line.rstrip().split("\t")
				new=line[5:len(line)]
				a=[i for i,j in enumerate(new) if j!="R"]
				#b=[j for i,j in enumerate(new) if j!="R"] 
				if len(a)==1:
					s=a[0]+5
					ss=line[0:5]
					ss.append(line[s])
					if len(hd[s].split("|"))>1:
						o.write(hd[s].split("|")[0]+"\t"+hd[s].split("|")[1]+"\t"+"\t".join(ss)+"\n")
					else:
						o.write(hd[s].split("|")[0]+"\t"+"Unknown"+"\t"+"\t".join(ss)+"\n")

	f.close()
	o.close()

def custom_annotation(in_file,in_anno,out_dir,out_pref):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form2.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form2.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form2.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form2.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form2.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form2.txt")
	ln,hd,tabgmut,tabcg,tabpg,tabv,tabef,tab = read_in_table_form2(in_file)
	sanno={}
	with open(in_anno) as tmp:
		for line in tmp:
			line=line.rstrip().split("\t")
			for k in tab.keys():
				if line[0]==k[0]:
					if int(k[1])>=int(line[1]) and int(k[1])<=int(line[2]):
						if k in sanno.keys():
							sanno[k].append(line[3])
						else: 
							sanno[k]=[line[3]]
	for k in sanno.keys():
		sanno[k]=list(set(sanno[k]))
	for k in tab.keys():
		if k in sanno.keys():
			tab[k]=["|".join(sanno[k])]+tab[k]
		else:
			tab[k]=["NULL"]+tab[k]
	for k in tabgmut.keys():
		if k in sanno.keys():
			tabgmut[k]=["|".join(sanno[k])]+tabgmut[k]
		else:
			tabgmut[k]=["NULL"]+tabgmut[k]
	for k in tabcg.keys():
		if k in sanno.keys():
			tabcg[k]=["|".join(sanno[k])]+tabcg[k]
		else:
			tabcg[k]=["NULL"]+tabcg[k]
	for k in tabpg.keys():
		if k in sanno.keys():
			tabpg[k]=["|".join(sanno[k])]+tabpg[k]
		else:
			tabpg[k]=["NULL"]+tabpg[k]
	for k in tabv.keys():
		if k in sanno.keys():
			tabv[k]=["|".join(sanno[k])]+tabv[k]
		else:
			tabv[k]=["NULL"]+tabv[k]
	for k in tabef.keys():
		if k in sanno.keys():
			tabef[k]=["|".join(sanno[k])]+tabef[k]
		else:
			tabef[k]=["NULL"]+tabef[k]
	hd2=["Custom_Anno"]+hd
	write_tabs_form2(tabcg,tabcg_out_table,hd2)
	write_tabs_form2(tabpg,tabpg_out_table,hd2)
	write_tabs_form2(tabv,tabv_out_table,hd2)
	write_tabs_form2(tabef,tabef_out_table,hd2)
	write_tabs_form2(tab,tab_out_table,hd2)
	write_tabs_form2(tabgmut,tabgmut_out_table,hd2)

def mask_seqs(in_file,in_msk,out_dir,out_pref):
	tabcg_out_table=os.path.join(out_dir,out_pref+".cds_change.form2.txt")
	tabpg_out_table=os.path.join(out_dir,out_pref+".prot_change.form2.txt")
	tabv_out_table=os.path.join(out_dir,out_pref+".mut_type.form2.txt")
	tabef_out_table=os.path.join(out_dir,out_pref+".effect.form2.txt")
	tab_out_table=os.path.join(out_dir,out_pref+".all.form2.txt")
	tabgmut_out_table=os.path.join(out_dir,out_pref+".gene_mut.form2.txt")
	ln,hd,tabgmut,tabcg,tabpg,tabv,tabef,tab = read_in_table_form2(in_file)
	idx=[]
	with open(in_msk) as tmp:
		for line in tmp:
			line=line.rstrip().split("\t")
			for k in tab.keys():
				if line[0]==k[0]:
					if int(k[1])>=int(line[1]) and int(k[1])<=int(line[2]):
						idx.append(k)
	idx=list(set(idx))
	tabn = dict([(key, val) for key, val in tab.items() if key not in idx]) 
	tabcgn = dict([(key, val) for key, val in tabcg.items() if key not in idx]) 
	tabpgn = dict([(key, val) for key, val in tabpg.items() if key not in idx]) 
	tabvn = dict([(key, val) for key, val in tabv.items() if key not in idx]) 
	tabefn = dict([(key, val) for key, val in tabef.items() if key not in idx]) 
	tabgmutn = dict([(key, val) for key, val in tabgmut.items() if key not in idx]) 
	write_tabs_form2(tabcgn,tabcg_out_table,hd)
	write_tabs_form2(tabpgn,tabpg_out_table,hd)
	write_tabs_form2(tabvn,tabv_out_table,hd)
	write_tabs_form2(tabefn,tabef_out_table,hd)
	write_tabs_form2(tabn,tab_out_table,hd)
	write_tabs_form2(tabgmutn,tabgmut_out_table,hd)

def mutation_summary(out_log,fl,fid,out,fm):
	with open(out,'w') as o:
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			a=[i for i,j in enumerate(hd) if fid in j]
			if len(a)==0:
				scn_print("Error: The selected ID is not in the summary table!")
				log_print(out_log,"Error: The selected ID is not in the summary table!")
				exit(1)
			elif len(a)>1:
				scn_print("Warning: This ID is duplicated in the summary table! Only the first one will be used here.")
				log_print(out_log,"Warning: This ID is duplicated in the summary table! Only the first one will be used here.")
				s=a[0]
			else:
				s=a[0]
			hd2=hd[4:s]+hd[(s+1):]
			if fm=="l":
				o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tmutation\tmutation_type\t"+"\t".join(hd2)+"\n")
			else:
				o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tmutation\tmutation_type\n")
			for line in f:
				line = line.rstrip().split("\t")
				#sam=line[s]
				if line[s]!="R":
					muts=line[4:s]+line[(s+1):]
					b=[i for i,j in enumerate(muts) if j!="R"]
					if len(b)==0:
						if fm=="l":
							if len(hd[s].split("|"))>1:
								o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Unique_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Unique_mutation"+"\t"+"\t".join(muts)+"\n")
						else:
							if len(hd[s].split("|"))>1:
								o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Unique_mutation\n")
							else:
								o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Unique_mutation\n")
					else:
						c=[i for i,j in enumerate(muts) if j==line[s]]
						if len(c)==0:
							if fm=="l":
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"New_mutation"+"\t"+"\t".join(muts)+"\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"New_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"New_mutation\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"New_mutation\n")
						else:
							if fm=="l":
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Known_mutation"+"\t"+"\t".join(muts)+"\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Known_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Known_mutation\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[s]+"\t"+"Known_mutation\n")
	f.close()
	o.close()

def mutation_summary_wanno(out_log,fl,fid,out,fm):
	with open(out,'w') as o:
		with open(fl,'r') as f:
			hd=f.readline()
			hd=hd.rstrip().split("\t")
			a=[i for i,j in enumerate(hd) if fid in j]
			if len(a)==0:
				scn_print("Error: The selected ID is not in the summary table!")
				log_print(out_log,"Error: The selected ID is not in the summary table!")
				exit(1)
			elif len(a)>1:
				scn_print("Warning: This ID is duplicated in the summary table! Only the first one will be used here.")
				log_print(out_log,"Warning: This ID is duplicated in the summary table! Only the first one will be used here.")
				s=a[0]
			else:
				s=a[0]
			hd2=hd[5:s]+hd[(s+1):]
			if fm=="l":
				o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tcustom_annot\tmutation\tmutation_type\t"+"\t".join(hd2)+"\n")
			else:
				o.write("ID\tregion\tchr\tpos\tgene_id\tgene_name\tcustom_annot\tmutation\tmutation_type\n")
			for line in f:
				line = line.rstrip().split("\t")
				#sam=line[s]
				if line[s]!="R":
					muts=line[5:s]+line[(s+1):]
					b=[i for i,j in enumerate(muts) if j!="R"]
					if len(b)==0:
						if fm=="l":
							if len(hd[s].split("|"))>1:
								o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Unique_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Unique_mutation"+"\t"+"\t".join(muts)+"\n")
						else:
							if len(hd[s].split("|"))>1:
								o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Unique_mutation\n")
							else:
								o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Unique_mutation\n")
					else:
						c=[i for i,j in enumerate(muts) if j==line[s]]
						if len(c)==0:
							if fm=="l":
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"New_mutation"+"\t"+"\t".join(muts)+"\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"New_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"New_mutation\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"New_mutation\n")
						else:
							if fm=="l":
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Known_mutation"+"\t"+"\t".join(muts)+"\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Known_mutation"+"\t"+"\t".join(muts)+"\n")
							else:
								if len(hd[s].split("|"))>1:
									o.write(fid+"\t"+hd[s].split("|")[1]+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Known_mutation\n")
								else:
									o.write(fid+"\t"+"UnKnown"+"\t"+line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[s]+"\t"+"Known_mutation\n")
	f.close()
	o.close()

def calculate_frequencies(fun_out_dir,fun_out_pref,fun_type):
	import pandas as pd
	in_table=os.path.join(fun_out_dir,fun_out_pref+"."+fun_type+".form2.txt")
	out_s=os.path.join(fun_out_dir,fun_out_pref+"."+fun_type+".sequence.frequencies.txt")
	out_p=os.path.join(fun_out_dir,fun_out_pref+"."+fun_type+".pos.frequencies.txt")
	fr=None
	fr2=None
	df = pd.read_csv(in_table, sep = '\t', \
	comment = '#', low_memory = False, header = 0)
	hd=list(df.columns)
	for i in hd[4:]:
		if fr is None:
			fr=df[i].value_counts().rename_axis('Mutation').reset_index(name='Count')
			fr.insert(0, "ID", i, True)
			fr.to_csv(out_s, index=False)
		else:
			nd=df[i].value_counts().rename_axis('Mutation').reset_index(name='Counts')
			nd.insert(0, "ID", i, True) 
			nd.to_csv(out_s, mode='a', header=False, index=False)
	for j in range(df.shape[0]):
		if fr2 is None:
			fr2=df.iloc[j][4:].value_counts().rename_axis('Mutation').reset_index(name='Counts')
			fr2.insert(0,"Chr",df.iloc[j][0],True)
			fr2.insert(1,"Pos",df.iloc[j][1],True)
			fr2.insert(2,"Gid",df.iloc[j][2],True)
			fr2.insert(3,"Gname",df.iloc[j][3],True)
			fr2.to_csv(out_p, index=False)
		else:
			nd=df.iloc[j][4:].value_counts().rename_axis('Mutation').reset_index(name='Counts')
			nd.insert(0,"Chr",df.iloc[j][0],True)
			nd.insert(1,"Pos",df.iloc[j][1],True)
			nd.insert(2,"Gid",df.iloc[j][2],True)
			nd.insert(3,"Gname",df.iloc[j][3],True)
			nd.to_csv(out_p, mode='a', header=False, index=False)

def summary_table_calculate_frequencies(out_dir,out_pref):
	import pandas as pd
	calculate_frequencies(out_dir,out_pref,"cds_change")
	calculate_frequencies(out_dir,out_pref,"prot_change")
	calculate_frequencies(out_dir,out_pref,"mut_type")
	calculate_frequencies(out_dir,out_pref,"effect")
	calculate_frequencies(out_dir,out_pref,"all")
	calculate_frequencies(out_dir,out_pref,"gene_mut")