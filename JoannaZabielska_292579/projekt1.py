#!/usr/bin/python
# -*- coding: utf-8 -*-
from Bio import SeqIO
import subprocess as sub

#INPUT
input_path = "./input/Input"
input_path2 = "./output/not_VC-mRNA.fas"
input_path3 = "./output/not_BC-mRNA.fas"
input_path4 = "./output/not_RA-mRNA.fas"
input_path5 = "./output/not_AA-mRNA.fas"

#DATABASES
univec_path = "./bazy/UniVec"
bacteria_path = "./bazy/BacteriaDB"
ref_path = "./bazy/RefProtDB"
all_prot_path = "./bazy/All_Protein_Database.fasta"
CDD_path = "./bazy/Cdd_LE/Cdd"

#MUSEQBOX
museqbox_path = "./MUSEQBOX/bin/MuSeqBox"



#FUNCTION

def seqnames_from_museqboxfile(filename):
  """wczytywanie identyfikatorow sekwencji z output'u programu MuSeqBox"""
  
  plik_msb = open(filename,"r")

  j=0
  nazwy_seq = []
  for i in plik_msb.readlines():

    if i[0]=="-" and j==0: 
	j+=1
	continue
    if i[0]=="-" and j>0: 
      break
    
    if j>0 and i[0] != "\n" and i[0] != " ": 
      nazwy_seq.append( i.split()[0] )

  return nazwy_seq

  
def MuSeqBox_Partition(museqboxfile_path, inputfile_path):

  museqboxfile = open(museqboxfile_path,"r")
  
  l = seqnames_from_museqboxfile(museqboxfile_path)
  
  sth_mRNA = []
  not_sth_mRNA = []
  
  for record in SeqIO.parse(open(inputfile_path, "rU"), "fasta"):

    if record.name.split("|")[-1] in l:
      sth_mRNA.append(record)
    else:
      not_sth_mRNA.append(record)

  return (sth_mRNA, not_sth_mRNA)

#MUSEQBOX to step 5
def blast_5(plik_blast, plik, kol, dobre,zle):
	import sys
	
	kol = int(kol)-1
	assert kol >= 0
	
	ids = set()
	x =open(plik_blast, "rU")
	for line in x:
		ids.add(line.split("\t")[kol])
	x.close()
	
	wczytane=open(plik,"rU")
	dobre_wczytane = open(dobre,"w")
	zle_wczytane = open(zle,"w")

	m,n=[],[]

	for record in SeqIO.parse(wczytane,"fasta"):
		if record.name in ids:
			m.append(record)
			SeqIO.write(record,dobre_wczytane,"fasta")
		else:
			n.append(record)
			SeqIO.write(record,zle_wczytane,"fasta")
	wczytane.close()
	dobre_wczytane.close()
	zle_wczytane.close()
	return(m,n)
				


#STEP 1

komenda1_1 = "makeblastdb -dbtype nucl -in %s -parse_seqids" % (univec_path)
komenda1_2 = "blastn -task blastn -db %s -query %s -evalue 1e-25 -show_gis -best_hit_overhang 0.25 -best_hit_score_edge 0.05 -out ./output/blastn_Vector" % (univec_path, input_path)
komenda1_3 = "%s -i ./output/blastn_Vector > ./output/VC.msb" %(museqbox_path)

p = sub.Popen(komenda1_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

p1 = sub.Popen(komenda1_2,shell=True, stdout=sub.PIPE)
p1.stdout.readlines()

p2 = sub.Popen(komenda1_3,shell=True, stdout=sub.PIPE)
p2.stdout.readlines()

(VC_mRNA, not_VC_mRNA) =  MuSeqBox_Partition("./output/VC.msb", input_path)
len_VC_mRNA =  len(VC_mRNA)

output_handle1 = open("./output/VC-mRNA.fas", "w")
SeqIO.write(VC_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_VC-mRNA.fas", "w")
SeqIO.write(not_VC_mRNA, output_handle2, "fasta")
output_handle2.close()
    

#STEP 2

komenda2_1 = "makeblastdb -dbtype nucl -in %s -parse_seqids" % (bacteria_path)
komenda2_2 = "blastn -task blastn -db %s -query %s -evalue 1e-30 -show_gis -best_hit_overhang 0.1 -best_hit_score_edge 0.05  -out ./output/blastn_Bacteria" % (bacteria_path, input_path2)
komenda2_3 = "%s -i ./output/blastn_Bacteria > ./output/BC.msb" %(museqbox_path)

p = sub.Popen(komenda2_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

p1 = sub.Popen(komenda2_2,shell=True, stdout=sub.PIPE)
p1.stdout.readlines()

p2 = sub.Popen(komenda2_3,shell=True, stdout=sub.PIPE)
p2.stdout.readlines()

(BC_mRNA, not_BC_mRNA) =  MuSeqBox_Partition("./output/BC.msb", input_path2)
len_BC_mRNA =  len(BC_mRNA)

output_handle1 = open("./output/BC-mRNA.fas", "w")
SeqIO.write(BC_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_BC-mRNA.fas", "w")
SeqIO.write(not_BC_mRNA, output_handle2, "fasta")
output_handle2.close()


#STEP 3

komenda3_1 = "makeblastdb -dbtype prot -in %s -parse_seqids" % (ref_path)
komenda3_2 = "blastx -db %s -query %s -evalue 1e-20 -show_gis -best_hit_overhang 0.1 -best_hit_score_edge 0.01 -out ./output/blastx_RefProt" % (ref_path, input_path3)
komenda3_3 = "%s -i ./output/blastx_RefProt > ./output/RA.msb" %(museqbox_path)

p = sub.Popen(komenda3_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

p1 = sub.Popen(komenda3_2,shell=True, stdout=sub.PIPE)
p1.stdout.readlines()

p2 = sub.Popen(komenda3_3,shell=True, stdout=sub.PIPE)
p2.stdout.readlines()

(RA_mRNA, not_RA_mRNA) =  MuSeqBox_Partition("./output/RA.msb",input_path3)
len_RA_mRNA = len(RA_mRNA)

output_handle1 = open("./output/RA-mRNA.fas", "w")
SeqIO.write(RA_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_RA-mRNA.fas", "w")
SeqIO.write(not_RA_mRNA, output_handle2, "fasta")
output_handle2.close()

#STEP 3.1

komenda3_1_1 = "%s -i ./output/RA.msb  -n 1  -p 4 -s 1 -F 10 10 0 0 95.0 62.0 ./output/fullcds.msb" % (museqbox_path)

p = sub.Popen(komenda3_1_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

(FL_mRNA, not_FL_mRNA) =  MuSeqBox_Partition("./output/fullcds.msb","./output/RA-mRNA.fas")
len_fullcds = len(FL_mRNA)
len_not_fullcds = len(not_FL_mRNA)

output_handle1 = open("./output/FL-mRNA.fas", "w")
SeqIO.write(FL_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_FL-mRNA.fas", "w")
SeqIO.write(not_FL_mRNA, output_handle2, "fasta")
output_handle2.close()

#STEP 3.2

komenda3_2_1 = "%s -i ./output/RA.msb -n 1  -p 4 -s 1 -F 10 10 0 0 60.0 62.0 ./output/PC.msb" % (museqbox_path)

p = sub.Popen(komenda3_2_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

(PC_mRNA, not_PC_mRNA) =  MuSeqBox_Partition("./output/PC.msb","./output/not_FL-mRNA.fas")
len_PC_mRNA = len(PC_mRNA)
len_not_PC_mRNA = len(not_PC_mRNA)

output_handle1 = open("./output/PC_mRNA.fas", "w")
SeqIO.write(PC_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_PC_mRNA.fas", "w")
SeqIO.write(not_PC_mRNA, output_handle2, "fasta")
output_handle2.close()


#STEP 4

komenda4_1 = "makeblastdb -dbtype prot -in %s -parse_seqids" % (all_prot_path)     
komenda4_2 = "blastx -db %s -query %s -evalue 1e-30 -show_gis -best_hit_overhang 0.1 -best_hit_score_edge 0.05 -out ./output/blastx_AllProt" % (all_prot_path, input_path4)
komenda4_3 = "%s -i ./output/blastx_AllProt > ./output/AA.msb" %(museqbox_path)

p = sub.Popen(komenda4_1,shell=True, stdout=sub.PIPE)
p.stdout.readlines()

p1 = sub.Popen(komenda4_2,shell=True, stdout=sub.PIPE)
p1.stdout.readlines()

p2 = sub.Popen(komenda4_3,shell=True, stdout=sub.PIPE)
p2.stdout.readlines()

(AA_mRNA, not_AA_mRNA) =  MuSeqBox_Partition("./output/AA.msb",input_path4)
len_AA_mRNA =  len(AA_mRNA)

output_handle1 = open("./output/AA-mRNA.fas", "w")
SeqIO.write(AA_mRNA, output_handle1, "fasta")
output_handle1.close()

output_handle2 = open("./output/not_AA-mRNA.fas", "w")
SeqIO.write(not_AA_mRNA, output_handle2, "fasta")
output_handle2.close()

#STEP 5

    
komenda5_2 = "rpstblastn -db %s -query %s  -evalue 1e-10  -out ./output/rpstblastn_ProtDomain -outfmt 6" % (CDD_path, input_path5)

p1 = sub.Popen(komenda5_2,shell=True, stdout=sub.PIPE)
p1.stdout.readlines()

CD_mRNA, not_CD_mRNA =  blast_5("./output/rpstblastn_ProtDomain", input_path5, 1, "./output/CD-mRNA.fas","./output/not_CD-mRNA.fas")

len_CD_mRNA = len(CD_mRNA)
len_not_CD_mRNA = len(not_CD_mRNA)

#STEP 6


inp = len(list(SeqIO.parse(open(input_path, "rU"), "fasta")))

  
summary_text = """

mRNAmarkup Report

  Number of input sequences:                            %d
  Number of potential vector-contaminated sequences:    %d (file: VC_mRNA.fas)
  Number of potential bacterial-contaminated sequences: %d (file: BC_mRNA.fas)
  Number of sequences matching the ReferenceDB:         %d (file: RA_mRNA.fas)
    Number of potential full-length coding sequences:     %d (file: FL_mRNA.fas)
    Non-qualifying sequences:                             %d (file: not_FL_mRNA.fas)
      Number of potential chimeric sequences:               %d (file: PC_mRNA.fas)
      Non-qualifying sequences:                             %d (file: not_PC_mRNA.fas)
  Number of sequences matching the AllProteinDB:        %d (file: AA_mRNA.fas)
  Number of sequences matching the ProteinDomainDB:     %d (file: CD_mRNA.fas)
  Number of remaining sequences:                        %d (file: remaining-mRNA)
 
""" % (inp, len_VC_mRNA, len_BC_mRNA, len_RA_mRNA, len_fullcds, len_not_fullcds, len_PC_mRNA, len_not_PC_mRNA, len_AA_mRNA, len_CD_mRNA,inp - (len_VC_mRNA + len_BC_mRNA + len_RA_mRNA  + len_AA_mRNA + len_CD_mRNA) )


summary = open("summary.txt","w")
summary.write(summary_text)
summary.close()

