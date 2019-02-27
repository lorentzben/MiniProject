Comp Bio 383 Mini Project Python Wrapper Option 1
-------------------------------------------------
E. coli strain K-12 was the first strain to have its genome sequenced. However it is not the only strain that researchers use.
As time goes on, the strains used by researches will evolve and as such some researchers have begun resequencing projects to see how the strain has evolved from the original sample. 
This program aims to automate the collection, assembly, and annotation of the E. coli K-12 sequences. 

## Prerequisities
* Python 3
* Linux
* SPAdes v3.11.1
* bowtie v1.2.2
* tophat v2.1.1
* cufflinks v2.2.1 
## Install

```shell
$ git clone https://github.com/lorentzben/MiniProject.git
```

After cloning a folder called MiniProject will be created. Inside will be two .py files and this README. 

## Running 
```shell
$ python3 main.py
```
This script will run automatically in a directory called Ben_Lorentz_Option_A. The data nessecary for analysis will be downloaded from NCBI's databases. 
Running
```shell
$ python3 cleanup.py
```
will delete the directory created by main.py. This is useful in debugging and when modifying code to be uploaded to github. 
## Output
Inside of the directory Ben_Lorentz_Option_A there will be OptionA.log which will contain parameters used as well as pertinant information. 
The second file of interest is Option1.fpkm which contains the output from cufflinks formatting as a .csv

## Version
* Version 1.0

## Author
* Ben Lorentz

## Future Plans
* refactor variable names to be more pythonic
* refactor to each chunk written in functions
* evaluate the use of subprocess as opposed to os.system, in an effort to pip st.out for verbose logging
* generalize input using accession numbers 
