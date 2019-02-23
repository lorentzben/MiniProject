#Automate Ecoli
#Ben Lorentz 2.22.19
import sys
import os
import logging
dir = os.getcwd()
cwd = dir+"/OptionA_Ben_Lorentz"
if(not os.path.exists("OptionA_Ben_Lorentz")):
    os.system("mkdir OptionA_Ben_Lorentz")
    os.chdir(dir+"/OptionA_Ben_Lorentz")
else:
    os.chdir(dir+"/OptionA_Ben_Lorentz")
resultdict = {}
longBois = []
logging.basicConfig(level=logging.DEBUG, filename="OptionA.log")
sraAccess = "SRR8185310"
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
if(not os.path.isfile("longBoiContigs.fa")):
    result = ""
    for item in longBois:
        result = result+ ">" +item[0].strip()+ '\n'
        result = result + item[1].strip() + '\n'
    #Boilerplate to write to file
    file = open("longBoi.fasta", "w")
    file.write(str(result))
    file.close()
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
print("ok")
