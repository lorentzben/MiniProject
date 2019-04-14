#!/usr/bin/python3 
"""Main.py: This is a script that automates the collection, assembly, and annotation of the E. coli K-12 sequences."""

__author__ = "Ben Lorentz"
__email__ = "blorentz@luc.edu"

import subprocess
import sys
import os
import logging


#runs bash command but does not print out
def normal(cmd):
    subprocess.run([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
    logging.info(cmd)
#runs call and prints out to shell session
def normal_verbose(cmd):
    call = subprocess.run([cmd], check=True, shell=True)
    logging.info(cmd)

#checks to see if working dir exists then moves in, or creates it based on necessity
def create_working_directory():
    if(not os.path.exists("OptionA_Ben_Lorentz")):
        os.system("mkdir OptionA_Ben_Lorentz")
        os.chdir(dir+"/OptionA_Ben_Lorentz")
    else:
        os.chdir(dir+"/OptionA_Ben_Lorentz")

#will download the .sra files if not found in the directory
def collect_sra_files():
    print("Getting Files")
    logging.info("The reads were not found so we are downloading them now")
    normal("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR818/SRR8185310/SRR8185310.sra")

#Converts .sra files into .fastq files
def sra_to_fastq():
    print("Turning .sra to .fastq")
    logging.info("Turning .sra file into .fastq file")
    normal("fastq-dump -I --split-files "+sraAccess+".sra")

#Assembles the reads using typical params with spades
def spades():
    print("Running spades with standard params")
    logging.info("spades -k 55,77,99,127 -t 2 --only-assembler -s "+sraAccess+"_1.fastq -o "+os.getcwd()+"/spades")
    os.system("spades -k 55,77,99,127 -t 2 --only-assembler -s "+sraAccess+"_1.fastq -o "+os.getcwd()+"/spades") 

#Goes into  spades folder and pull contigs.fasta out
if (not os.path.isfile("contigs.fasta")):
    os.chdir(cwd+"/spades")
    os.system("cp contigs.fasta "+cwd+"/")
    os.chdir(cwd)
if(resultdict == {}):
    with open("contigs.fasta", "r") as myfile:
        allseq = myfile.readlines()
    #My way of parsing a fasta file into a dictionary keys are seqID and values are seqs
    seq = ""
    header = ""
    for item in allseq:
        if(item[0]==">"):
             header =  item[1:].strip()
             seq = " "
        if (item[0] != ">" and (item[0] == "A" or item[0] == "T" or item[0] == "C" or item[0] == "G")):
             seq = seq + item[0:].strip()
        resultdict[header]=seq

    #this turns a dict into a list
    resultlist = list(resultdict.items())
    logging.info(str(len(resultlist)) + " contigs after alignment")
    for item in resultlist:
        if (len(item[1]) >= 1000):
            longBois.append(item)
    logging.info("There are " + str(len(longBois)) + " contigs > 1000 in the assembly.")
#Calculate the length of the assembly 
assemb = 0
if(assemb == 0):
    for each in longBois:
        assemb = len(each[1]) + assemb
    logging.info("There are " +str(assemb)+ " bp in the assembly")
#Pulls in the contigs over 1000 bp
if(not os.path.isfile("longBoiContigs.fa")):
    result = ""
    for item in longBois:
        result = result+ ">" +item[0].strip()+ '\n'
        result = result + item[1].strip() + '\n'
    #Boilerplate to write to file
    file = open("longBoi.fasta", "w")
    file.write(str(result))
    file.close()
#Runs prokka on the contigs longer than 1000 bp 
if(not os.path.exists(os.getcwd()+"/prokka")):
    os.system("prokka --usegenus --outdir prokka -genus Escherichia longBoi.fasta")
    logging.info("prokka --usegenus --outdir prokka -genus Escherichia longBoi.fasta")
    os.chdir(cwd+"/prokka")
    os.system("cp PROKKA*.txt "+cwd+"/prokka.txt")
    os.chdir(cwd)
    with open("prokka.txt", "r") as prokka:
        allLine = prokka.readlines()
    for line in allLine:
        logging.info(line.strip())
ecoli = {}
#Calculating diffs from the reference sequence and writes to log
if(ecoli == {}):
    ecoli = {"CDS":4140, "tRNAs":89}
    with open("prokka.txt", "r") as prokka:
       allLine = prokka.readlines()
    data = {"CDS":int(allLine[3].split(":")[1].strip()), "tRNAs": int(allLine[4].split(":")[1].strip())}
    cdsDiff=ecoli["CDS"]-data["CDS"]
    tRNADiff=ecoli["tRNAs"]-data["tRNAs"]
    logging.info("Prokka found " +str(abs(cdsDiff))+ (" additional" if cdsDiff > 0 else " fewer") + " CDS and " +str(abs(tRNADiff))+ (" additional " if  tRNADiff > 0 else " fewer ") + "tRNA than the RefSeq.")  
#Pulls down sequence of ecoli k-12 and makes index and calls tophat on the reference
if(not os.path.isfile("SRR1411276.sra")):
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR141/SRR1411276/SRR1411276.sra")
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna")
    print("Turning .sra to .fastq")
    logging.info("Turning .sra file into .fastq file")
    os.system("fastq-dump -I --split-files SRR1411276.sra")
    logging.info("fastq-dump -I --split-files SRR1411276.sra")
    os.system("bowtie2-build NC_000913.fna EcoliK12")
    logging.info("bowtie2-build NC_000913.fna EcoliK12")
    os.system("tophat2 --no-novel-juncs -o tophat EcoliK12 SRR1411276_1.fastq")
    logging.info("tophat2 --no-novel-juncs -o tophat EcoliK12 SRR1411276_1.fastq")
os.chdir("tophat")
#runs cufflinks on the output from tophat
if( not os.path.exists("cufflinks")):
    logging.info("running cufflinks")
    logging.info("cufflinks accepted_hits.bam -o cufflinks")
    os.system("cufflinks accepted_hits.bam -o cufflinks")
    os.chdir("cufflinks")
    os.system("cp transcripts.gtf" + cwd)
    os.chdir(cwd)
#moves the cufflinks output to the top folder
if( not os.path.isfile("transcripts.gtf")):
    logging.info("Moving the transcripts.gtf file to: " + cwd)
    os.chdir(cwd)
    os.chdir("tophat")
    os.chdir("cufflinks")
    os.system("cp transcripts.gtf " +  cwd)
    os.chdir(cwd)
#creates the Option1.fpkm file formatted correctly
if( not os.path.isfile("Option1.fpkm")):
    logging.info("Formatting the final output")
    output = ""
    with open("transcripts.gtf", "r") as transcript:
        record = []
        for line in transcript:
            record.append(line)
        for i in range(0,len(record)):
            if(i % 2 == 0):
                seqname = record[i].split('\t')[0].strip()
                start = record[i].split('\t')[3].strip()
                end = record[i].split('\t')[4].strip()
                strand = record[i].split('\t')[6].strip()
                FPKM = record[i].split('\t')[8].split(";")[2]
            else:
                seqname = record[i].split('\t')[0].strip()
                start = record[i].split('\t')[3].strip()
                end = record[i].split('\t')[4].strip()
                strand = record[i].split('\t')[6].strip()
                FPKM = record[i].split('\t')[8].split(";")[3]
            output = output + seqname+ ',' +start+ "," +end+ ',' +strand+ "," +FPKM + '\n'
    
    #Boilerplate to write to file
    file = open("Option1.fpkm", "w")
    file.write(str(output))
    file.close()


def main():
    #TODO These need to go after all of the functions
    dir = os.getcwd()
    cwd = dir+"/OptionA_Ben_Lorentz"
    create_working_directory()
    #TODO These need to go after all of the functions
    result_dict = {}
    long_Bois = []
    logging.basicConfig(level=logging.DEBUG, filename="OptionA.log")
    sraAccess = "SRR8185310"
    if(not os.path.isfile(str(sraAccess+".sra"))):
       collect_sra_files()
    if(not os.path.isfile(str(sraAccess+"_1.fastq"))):
        sra_to_fastq() 
    if(not os.path.exists(os.getcwd()+"/spades")):
        spades()
    

print("The final file is called Option1.fpkm in the directory : " + cwd)
if __init__ == "__main__"":
    main()
