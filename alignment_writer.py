#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nom du projet : SeqRef
Script : alignment_writer.py

Description :
Ce script permet de générer une séquence de référence d'un transcrit selon son identifiant RefSeq.
Il fait usage de l'API NCBI REST API Datasets v2 pour récuperer les metadatas des transcrits et de l'API REST d'UCSC pour récupérer la séquence génomique.

Auteur : Romain LEVERGEOIS
Email : romain.levergeois@aphp.fr
Organisation : Hopital Armand Trousseau APHP - Service de Génétique
Date : 06/09/2024
Version : 1.0.0

Licence :
Ce script est publié sous la licence MIT.
Voir le fichier LICENSE pour plus de détails.

Dépendances :
Ce projet utilise API NCBI REST API Datasets v2, développé par le National Institutes of Health/Department of Health and Human Services, Bethesda, MD, U.S.A. Tous droits réservés.
Ce projet utilise UCSC REST API, développé par le Genomics Institute de l'Université de Californie Santa Cruz (UCSC),  Tous droits réservés.
Veuillez consulter le fichier LICENSE pour plus de détails sur la redistribution et les conditions d'utilisation.
"""

# Import libraries ##################################################################
import sys
import os
import math
import datetime
import re
import json 
import requests
from zipfile import ZipFile


# Functions #########################################################################
def encrypt(string,char,length):
    output=char.join(string[i:i+length] for i in range(0,len(string),length))
    return(output)

def get_metadata(accession):
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/accession/{}?returned_content=PRODUCT".format(accession)
    response = requests.get(url)
    
    # Check if the response is successful
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        return None


def get_sequence_hg38(orientation, chr, start, end):
	tmp=(1 if orientation == "minus" else 0)
	url = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{};start={};end={};revComp={}".format(chr,start,end, tmp)
	response = requests.get(url)

	# Check if the response is successful
	if response.status_code == 200:
		data = response.json()
		return data
	else:
		return None

def get_prots(gene_id, filename):
	url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/download?filename={}".format(filename)
	data = {"gene_ids": [gene_id], "include_annotation_type": ["FASTA_PROTEIN"]}
	response = requests.post(url, json = data, headers={'content-type': 'application/json','accept': 'application/zip'})

	# Check if the response is successful
	if response.status_code == 200:
		with open(filename+'.zip', 'wb') as file:
			file.write(response.content)

# BEGIN #############################################################################
## Variables
if (len(sys.argv)!=4):
	print("Erreur : nombre d'arguments insufisant")
	exit(1)

accession = sys.argv[1] 
CINQprime = int(sys.argv[2]) 
TROISprime = int(sys.argv[3]) 

metadata_prod = get_metadata(accession)

if (metadata_prod=={}):
	print("Erreur: NM non reconnu")
	exit(2)

if (int(metadata_prod["total_count"])==0 or metadata_prod=={}):
	print("Erreur: NM non reconnu")
	exit(2)	

geneSymbol=metadata_prod["reports"][0]["product"]["symbol"]
filenamezip=geneSymbol+".zip"
filename=geneSymbol


get_prots(metadata_prod["reports"][0]["product"]["gene_id"], filename)

## Read proteins zip file
with ZipFile(filenamezip) as myzip:
    with myzip.open("ncbi_dataset/data/protein.faa") as myfile:
        prots=str(myfile.read()).split('>')[1:]


## Reading all transcripts
for transcript in metadata_prod["reports"][0]["product"]["transcripts"]:
	if (str(transcript["accession_version"])!=accession):
		continue

	### Permet de gerer les chr alternatifs et recuperer que les exons pour GRCh38
	exons=0
	indexe=-1
	r = re.compile(".*GRCh38.p14 Primary Assembly.*")
	for i in range(0,len(transcript["genomic_locations"])):
		if re.match(".*GRCh38.p14 Primary Assembly.*",transcript["genomic_locations"][i]["sequence_name"]) :
			indexe=i

	if indexe!=-1 :
		genomic_locations=transcript["genomic_locations"][indexe]

		### Debut et Fin des exons
		exonStarts=[]
		exonEnds=[]
		for exon in genomic_locations["exons"]:
			exonStarts.append(int(exon["begin"]))
			exonEnds.append(int(exon["end"]))

		### Postions chromosomiques du gene
		start=genomic_locations["genomic_range"]["begin"]
		end=genomic_locations["genomic_range"]["end"]

		### Extraction chromosome du NC
		chr=int(genomic_locations["genomic_accession_version"].split('.')[0][-2:len(genomic_locations["genomic_accession_version"].split('.')[0])])
		if (chr==23):
			chr="X"
		if (chr==24):
			chr="Y"

		### Enregistre le NM
		name = transcript["accession_version"]


		## On recupere la sequence fasta avec la REST API UCSC 
		orientation=genomic_locations["genomic_range"]["orientation"]
		
		returned = get_sequence_hg38(orientation, chr, int(start)-CINQprime,int(end)+TROISprime)["dna"]

		## Mettre tous les carateres en minuscule
		fasta=returned.lower()
		#fasta=''.join(fasta.splitlines()[1:])
		#fasta=fasta.replace("\n","")


		## Position relative des exons
		tmpExonStarts=[]
		tmpExonEnds=[]

		if (orientation=="minus"):
			for eStart,eEnd in zip(exonStarts,exonEnds):
				tmpExonStarts.append(abs(eStart-int(end))+1)
				tmpExonEnds.append(abs(eEnd-int(end)))
				ExonEnds=tmpExonStarts
				ExonStarts=tmpExonEnds
		else:
			for eStart,eEnd in zip(exonStarts,exonEnds):
				tmpExonStarts.append(abs(eStart-int(start)))
				tmpExonEnds.append(abs(eEnd-int(start))+1)
				ExonEnds=tmpExonEnds
				ExonStarts=tmpExonStarts            


		### Correction of positions with +/- n bases in 3' and 5'.'
		newExonStarts=[]
		newExonEnds=[]
		for t in range(len(ExonStarts)):
			newExonStarts.append(ExonStarts[t]+CINQprime)

		for t in range(len(ExonEnds)):
			newExonEnds.append(ExonEnds[t]+CINQprime)

		## Construction of the fasta sequence with exons in uppercase
		deb=0
		fasta_exon=""
		fasta_exon=fasta_exon+fasta[0:newExonStarts[0]]
		for i in range(0,len(newExonStarts)-1):
			fasta_exon=fasta_exon+fasta[newExonStarts[i]:newExonEnds[i]].upper()+fasta[newExonEnds[i]:newExonStarts[i+1]]

		fasta_exon=fasta_exon+fasta[newExonStarts[len(newExonStarts)-1]:newExonEnds[len(newExonEnds)-1]].upper()+fasta[newExonEnds[len(newExonEnds)-1]:len(fasta)]

		## On converti la chaine de caracteres en tableau de caracteres
		tab_nucleotides=list(fasta_exon)

		type_nm=transcript["type"]
		if (type_nm!="NON_CODING"):
			### CDS
			cdsStart_rna=int(transcript["cds"]["range"][0]["begin"])
			cdsEnd_rna=int(transcript["cds"]["range"][0]["end"])

			## Relative positions of CDS
			### Begin
			tailleExon=0
			numExon_cdsStart=-1
			while(tailleExon<cdsStart_rna):
				numExon_cdsStart+=1
				tailleExon=tailleExon+ExonEnds[numExon_cdsStart]-ExonStarts[numExon_cdsStart]

			diff=tailleExon-(ExonEnds[numExon_cdsStart]-ExonStarts[numExon_cdsStart])
			pas=cdsStart_rna-diff
			cdsStart=ExonStarts[numExon_cdsStart]+pas
			cdsStart=cdsStart-1 # We wnat format M.. not ..M (-1 cause python  index bdegins at 0)
	

			### End
			tailleExon=0
			numExon_cdsEnd=-1
			while(tailleExon<cdsEnd_rna):
				numExon_cdsEnd+=1
				tailleExon=tailleExon+ExonEnds[numExon_cdsEnd]-ExonStarts[numExon_cdsEnd]

			diff=tailleExon-(ExonEnds[numExon_cdsEnd]-ExonStarts[numExon_cdsEnd])
			pas=cdsEnd_rna-diff
			cdsEnd=ExonStarts[numExon_cdsEnd]+pas
			cdsEnd=cdsEnd-1
		

			### Correction of positions with +/- n bases in 3' and 5'.'
			cdsStart=cdsStart+CINQprime
			cdsEnd=cdsEnd+CINQprime

			### Recupere la bonne sequence proteique
			NP=transcript["protein"]["accession_version"]

			for seq in prots:
				x = re.findall(NP, seq)
				if(x):
					prot=seq

			prot=''.join(prot.split('\\n')[1:])
			prot = prot[0:-1]
			prot=encrypt(prot,'..',1)+".."+"*"
			## Convert the string to amino acids table
			tab_aa=list(prot)

			## Synchronise les tableaux
			newTab_aa=[' ' for i in range(len(tab_nucleotides))]

			deb=0
			fin=0
			for o in range(numExon_cdsStart,numExon_cdsEnd+1):
				if (deb==0):
					taille=newExonEnds[o]-cdsStart
					fin=fin+taille
					for j,k in zip(range(cdsStart,newExonEnds[o]),tab_aa[deb:fin]):
						newTab_aa[j]= k
				else:
					taille=newExonEnds[o]-newExonStarts[o]
					fin=fin+taille
					for j,k in zip(range(newExonStarts[o],newExonEnds[o]),tab_aa[deb:fin]):
						newTab_aa[j]= k
				deb=fin
	else :
		print("Erreur: pas de reference hg38 pour ce transcrit")
		os.remove(filenamezip)
		exit(3)


### Format
deb=0
fin=60
final=""

for i in range(0,int(math.ceil(len(tab_nucleotides)/60))):

	final=final+(''.join(tab_nucleotides[deb:fin]))+'\n'
	if (type_nm!="NON_CODING"):
		final=final+(''.join(newTab_aa[deb:fin]))+'\n'
	deb=fin
	fin+=60
final=final+(''.join(tab_nucleotides[fin-60:len(tab_nucleotides)]))+'\n'

#### Header of aligment text

print("\t"+datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
print("\tGENE ETUDIE: "+geneSymbol+" - "+metadata_prod["reports"][0]["product"]["description"])
print("\tCHROMOSOME : "+str(chr))
print("\tBRIN : "+orientation)
print("\tTRANSCRIT : "+name)
print("\tNombre d'exons : "+str(len(exonStarts)))
print("\tPositions des exons :")
cpt=1
posExonString = ""
for a,b in zip(newExonStarts,newExonEnds):
	if(cpt==6):
		posExonString+=str(a+1)+"-"+str(b)+",\n"
		cpt=1
	elif(cpt==1):	
		posExonString+="\t\t"+str(a+1)+"-"+str(b)+", "
		cpt+=1
	else :
		posExonString+=str(a+1)+"-"+str(b)+", "
		cpt+=1
print(posExonString)
if (type_nm!="NON_CODING"):
	print("\tNum&eacute;ro des exons codants : "+str(numExon_cdsStart+1)+" a "+str(numExon_cdsEnd+1))
	print("\tPositions de la S&eacute;quence CoDante (CDS) : "+str(cdsStart+1)+" a "+str(cdsEnd-1)+"\n\n")

#### MAIN
deb_fasta=-60
fin_fasta=0
deb_prot=-19
fin_prot=0


for line in final.splitlines():
	if line.rstrip(' ').rstrip('\n'):
		x = re.findall("[ *.]", line)
		if (x):
			length=len(line.replace('.','').replace(' ',''))
			if length!=20:
				deb_prot=fin_prot+1
				fin_prot = deb_prot+length-1
				print(str(deb_prot)+"\t"+encrypt(line,' ',10)+"\t"+str(fin_prot)+"\n")
			else :
				deb_prot=fin_prot+1
				fin_prot=fin_prot+20
				print(str(deb_prot)+"\t"+encrypt(line,' ',10)+"\t"+str(fin_prot)+"\n")
		else :
			deb_fasta=fin_fasta+1
			fin_fasta=fin_fasta+len(line)
			print(str(deb_fasta)+"\t"+encrypt(line,' ',10)+"\t"+str(fin_fasta)+"\n")

os.remove(filenamezip)

# END #############################################################################
