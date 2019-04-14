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

def check_program(tool):
    try:
        subprocess.call([tool, "--help"],stdout=subprocess.PIPE)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            if tool == "fastq-dump":
                normal("sudo apt install -y sra-toolkit")
                logging.info("sudo apt install -y sra-toolkit")
            else:
                print(tool + " not found, installing")
                logging.info(tool + " not found")
                normal("sudo apt install -y " + tool)
                logging.info("sudo apt install -y " + tool)
    logging.info(tool +" installed and working fine")
        
def dependencies():
    depends = ["python3","spades","bowtie","tophat","cufflinks","fastq-dump"]
    for tool in depends:
        check_program(tool)

#checks to see if working dir exists then moves in, or creates it based on necessity
def create_working_directory(dir):
    if(not os.path.exists("OptionA_Ben_Lorentz")):
        os.system("mkdir OptionA_Ben_Lorentz")
        os.chdir(dir+"/OptionA_Ben_Lorentz")
    else:
        os.chdir(dir+"/OptionA_Ben_Lorentz")

#will download the .sra files if not found in the directory
def collect_sra_files():
    print("Getting Files")
    logging.info("The reads were not found so we are downloading them now")
    normal("fastq-dump -A SRR8185310 -O .")

#Converts .sra files into .fastq files
def sra_to_fastq():
    print("Turning .sra to .fastq")
    logging.info("Turning .sra file into .fastq file")
    normal("fastq-dump -I --split-files "+sraAccess+".sra")

#Assembles the reads using typical params with spades
def spades(sraAccess):
    print("Running spades with standard params")
    logging.info("spades -k 55,77,99,127 -t 2 --only-assembler -s "+sraAccess+".fastq -o "+os.getcwd()+"/spades")
    normal("spades -k 55,77,99,127 -t 2 --only-assembler -s "+sraAccess+".fastq -o "+os.getcwd()+"/spades") 

#Goes into  spades folder and pull contigs.fasta out
def pull_out_contigs(cwd):
    os.chdir(cwd+"/spades")
    normal("cp contigs.fasta "+cwd+"/")
    os.chdir(cwd)

#Calculates the number of contigs over 1000
def count_contigs(long_bois):
    with open("contigs.fasta", "r") as myfile:
        all_seq = myfile.readlines()
    #My way of parsing a fasta file into a dictionary keys are seqID and values are seqs
    seq = ""
    header = ""
    result_dict = {}
    for item in all_seq:
        if(item[0]==">"):
             header =  item[1:].strip()
             seq = " "
        if (item[0] != ">" and (item[0] == "A" or item[0] == "T" or item[0] == "C" or item[0] == "G")):
             seq = seq + item[0:].strip()
        result_dict[header]=seq

    #this turns a dict into a list
    result_list = list(result_dict.items())
    logging.info(str(len(result_list)).join(" contigs after alignment"))
    for item in result_list:
        if (len(item[1]) >= 1000):
            long_bois.append(item)
    logging.info("There are " + str(len(longBois)) + " contigs > 1000 in the assembly.")
    print(long_bois)
    return long_bois


#Calculate the length of the assembly 
def assembly_len(assemb, long_bois):
    for each in long_bois:
        assemb = len(each[1]) + assemb
    return assemb
    logging.info("There are " +str(assemb)+ " bp in the assembly") 

#Pulls in the contigs over 1000 bp
def write_fasta_to_file(long_bois):
    print("Working on writing longBoi.fasta")
    print(long_bois)
    result = ""
    for item in long_bois:
        result = result+ ">" +item[0].strip()+ '\n'
        result = result + item[1].strip() + '\n'
    #Boilerplate to write to file
    file = open("longBoi.fasta", "w")
    file.write(str(result))
    file.close()

#Runs prokka on the contigs longer than 1000 bp 
def prokka():
    normal("prokka --usegenus --outdir prokka -genus Escherichia longBoi.fasta")
    logging.info("prokka --usegenus --outdir prokka -genus Escherichia longBoi.fasta")
    os.chdir(cwd+"/prokka")
    normal("cp PROKKA*.txt "+cwd+"/prokka.txt")
    os.chdir(cwd)
    with open("prokka.txt", "r") as prokka:
        all_line = prokka.readlines()
    for line in all_line:
        logging.info(line.strip())

ecoli = {}
#Calculating diffs from the reference sequence and writes to log
def reference_diffs(ecoli):
    ecoli = {"CDS":4140, "tRNAs":89}
    with open("prokka.txt", "r") as prokka:
       all_line = prokka.readlines()
    data = {"CDS":int(all_line[3].split(":")[1].strip()), "tRNAs": int(all_line[4].split(":")[1].strip())}
    cds_diff=ecoli["CDS"]-data["CDS"]
    tRNA_diff=ecoli["tRNAs"]-data["tRNAs"]
    logging.info("Prokka found " +str(abs(cds_diff))+ (" additional" if cds_diff > 0 else " fewer") + " CDS and " +str(abs(tRNA_diff))+ (" additional " if  tRNA_diff > 0 else " fewer ") + "tRNA than the RefSeq.")  

#Pulls down sequence of ecoli k-12 and makes index and calls tophat on the reference
def ecoli_tophat():
    normal("fastq-dump -A SRR1411276 -O .")
    normal("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna")
    print("Turning .sra to .fastq")
    logging.info("Turning .sra file into .fastq file")
    normal("fastq-dump -I --split-files SRR1411276.sra")
    logging.info("fastq-dump -I --split-files SRR1411276.sra")
    normal("bowtie2-build NC_000913.fna EcoliK12")
    logging.info("bowtie2-build NC_000913.fna EcoliK12")
    normal("tophat2 --no-novel-juncs -o tophat EcoliK12 SRR1411276_1.fastq")
    logging.info("tophat2 --no-novel-juncs -o tophat EcoliK12 SRR1411276_1.fastq")


#runs cufflinks on the output from tophat
def cufflinks():
    logging.info("running cufflinks")
    logging.info("cufflinks accepted_hits.bam -o cufflinks")
    normal("cufflinks accepted_hits.bam -o cufflinks")
    os.chdir("cufflinks")
    normal("cp transcripts.gtf" + cwd)
    os.chdir(cwd)

#moves the cufflinks output to the top folder
def copy_cufflinks_out():
    logging.info("Moving the transcripts.gtf file to: " + cwd)
    os.chdir(cwd)
    os.chdir("tophat")
    os.chdir("cufflinks")
    normal("cp transcripts.gtf " +  cwd)
    os.chdir(cwd)

#creates the Option1.fpkm file formatted correctly
def write_fpkm_file():
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
    dir = os.getcwd()
    cwd = str(dir)+"/OptionA_Ben_Lorentz"

    create_working_directory(dir)

    long_bois = []
    logging.basicConfig(level=logging.DEBUG, filename="OptionA.log")
    sraAccess = "SRR8185310"

    dependencies()

    if(not os.path.isfile(str(sraAccess+".sra"))):
       collect_sra_files()
    if(not os.path.exists(os.getcwd()+"/spades")):
        spades(sraAccess)
    if (not os.path.isfile("contigs.fasta")):
        pull_out_contigs(cwd)
    if(long_bois == {}):
        long_bois  = count_contigs(long_bois)
        print(long_bois)
    assemb = 0 
    if(assemb == 0):
        assemb = assembly_len(assemb, long_bois)
    if(not os.path.isfile("longBoiContigs.fa")):
        print(long_bois)
        write_fasta_to_file(long_bois)
    if(not os.path.exists(os.getcwd()+"/prokka")):
        prokka() 
    ecoli = {}
    if(ecoli == {}):
        reference_diffs(ecoli)
    if(not os.path.isfile("SRR1411276.sra")):
        ecoli_tophat()
    os.chdir("tophat")
    if(not os.path.exists("cufflinks")):
        cufflinks()
    os.chdir(cwd)
    if(not os.path.isfile("transcripts.gtf")):
        copy_cufflinks_out() 
    if( not os.path.isfile("Option1.fpkm")):
        write_fpkm_out()
    print("The final file is called Option1.fpkm in the directory : " + cwd)

if __name__ == "__main__":
    main()

