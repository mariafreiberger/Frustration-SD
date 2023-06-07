import pandas as pd
import numpy as np
import os
import os.path
import matplotlib.pyplot as plt
from Bio import SeqIO
import sys

def frustration(pdb,chain):
	#------------ Genera los modelos y el msa -----------
	os.system('wget \'http://www.rcsb.org/pdb/files/'+pdb+'.pdb\' -O '+pdb+'.pdb')
	pdbfile=open(pdb+'.pdb')
	msa=open('msa.fasta','w')
	os.system('mkdir Models')
	directory=os.getcwd()+'/'
	n_models=0
	a=0
	for line in pdbfile.readlines():
		if line[0:5] == 'MODEL':
			sp=line.split()
			if a!=0:
				out.close()
				for record in SeqIO.parse(pdb+'.pdb', "pdb-atom"):
					if chain == record.annotations["chain"]:
						msa.write(">Model_"+sp[1]+"\n%s\n" % (record.seq))
			out=open('Models/Model_'+sp[1]+'.pdb','w')
			n_models=int(sp[1])
			a+=1
		if line[0:4] == 'ATOM' and line[21] == chain:
			out.write(line)
	out.close()
	pdbfile.close()
	for record in SeqIO.parse(pdb+'.pdb', "pdb-atom"):
		if chain == record.annotations["chain"]:
			msa.write(">Model_"+sp[1]+"\n%s\n" % (record.seq))
	msa.close()
	#------------ Calcula la frustración -------------
	rsc=open('Frusta.r','w')
	pdbdir=os.getcwd()+'/Models/'
	ResultsDir=os.getcwd()+'/'
	rsc.write('library(frustratometeR)\ndir_frustration(PdbsDir = \"'+pdbdir+'\", Mode = \"singleresidue\", ResultsDir = \"'+ResultsDir+'\")')
	rsc.close()
	os.system('Rscript Frusta.r')
	
	return n_models

def CreaDF(pdb,chain,models):
	a=0
	#-------------- DataFrames ------------
	for model in range(1,models+1):
		url='Model_'+str(model)+'.done/FrustrationData/Model_'+str(model)+'.pdb_singleresidue'
		df= pd.read_csv(url,sep=' ')
		id_col='fi_model_'+str(model)
		if a==0:
			df_total = pd.DataFrame(df['FrstIndex'])
			a+=1
		else:
			df_total.insert(1,id_col,df['FrstIndex'])
	return plots(df, df_total,pdb,chain)
	
def plots(df,df_total,pdb,chain):
	df_plot = pd.DataFrame({
		'res': df['Res'],
		'mean': df_total.mean(axis=1),
		'std': df_total.std(axis=1)})
	
	plt.errorbar(df_plot['res'],df_plot['mean'], df_plot['std'], linestyle='None', marker='o')
	plt.xlabel("Residue Number")
	plt.ylabel("Mean Frustration Index")
	plt.savefig(pdb+'_'+chain+'.png')
	plt.close()
	
#------------Main-------	
models=frustration(sys.argv[1],sys.argv[2])
CreaDF(sys.argv[1],sys.argv[2],models)
