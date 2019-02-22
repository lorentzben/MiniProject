#Automate Ecoli
#Ben Lorentz 2.22.19
import os
import logging
dir = os.getcwd()

if(not os.path.exists("OptionA_Ben_Lorentz")):
    os.system("mkdir OptionA_Ben_Lorentz")
    os.chdir(dir+"/OptionA_Ben_Lorentz")
else:
    os.chdir(dir+"/OptionA_Ben_Lorentz")
#Start Logging in the correct dir 
logging.basicConfig(filename='OptionA.log',level=logging.DEBUG)

#todo check that directory "depen" exists and has X Y and Z in it

if(not os.path.isfile("SRR8185310.sra")):
    logging.info("The reads were not found so we are downloading them now")
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR818/SRR8185310/SRR8185310.sra")

print("ok")
