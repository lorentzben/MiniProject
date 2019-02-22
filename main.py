#Automate Ecoli
#Ben Lorentz 2.22.19
import sys
import os
import logging
dir = os.getcwd()
if(not os.path.exists("OptionA_Ben_Lorentz")):
    os.system("mkdir OptionA_Ben_Lorentz")
    os.chdir(dir+"/OptionA_Ben_Lorentz")
else:
    os.chdir(dir+"/OptionA_Ben_Lorentz")

logging.basicConfig(level=logging.DEBUG, filename="OptionA.log")
sraAccess = "SRR8185310"
#TODO check that directory "depen" exists and has X Y and Z in it for PROKKA esp
#Downloads the files form sra database 
if(not os.path.isfile(str(sraAccess+".sra"))):
    print("Getting Files")
    logging.info("The reads were not found so we are downloading them now")
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR818/SRR8185310/SRR8185310.sra")
#Need to turn the .sra file into .fastq so we can use spades to assemble
if(not os.path.isfile(str(sraAccess+"_1.fastq"))):
    print("Turning .sra to .fastq")
    logging.info("Turning .sra file into .fastq file")
    os.system("fastq-dump -I --split-files "+sraAccess+".sra")

#Assemble reads using spades and typical params
if(not os.path.exists(os.getcwd()+"/spades")):
   print("Running spades with standard params")
   os.system("spades -k 55,77,99,127 -t 2 --only-assembler -s "+sraAccess+"_1.fastq -o "+os.getcwd()+"/spades")
#TODO need to go into the spades folder and pull contigs.fasta out
print("ok")
