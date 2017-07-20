#! /usr/bin/python
#
# Summary: retrieves a genome from the UCSC Genome Browser
#
# example call: 
# python downloadGenome.py --genome=mm9 --genomeDir=.

import sys
import ftplib
import os

def listGenomes():
    ftp = ftplib.FTP("hgdownload.cse.ucsc.edu", "ftp", "ftp")
    print("Available genomes: ")
    ftp.retrlines('LIST /goldenPath')
    

def downloadGenome(genome, genomeDir):
    os.chdir(genomeDir)    
    ftp = ftplib.FTP("hgdownload.cse.ucsc.edu", "ftp", "ftp")
    availFiles = []
    ftp.retrlines('LIST /goldenPath/'+genome+'/bigZips', availFiles.append)
    found = ""
    for availFile in availFiles:
        for file in ["chromFa.zip", "chromFa.tar.gz", genome+".zip", genome+".tar.gz", ""]:        
            if file != "" and availFile.find(file) >= 0: 
                found = file
    if found == "":
        print("No genome file found. Available files: "+str(availFiles))
    else:
        # download and unpack
        print("Downloading file: "+found)
        ftp.retrbinary("RETR /goldenPath/"+genome+"/bigZips/"+found, open(found, "wb").write)
        cmd = ""
        tokens = os.path.splitext(found)
        if tokens[len(tokens)-1] == "zip": cmd = "unzip "+found+" -d "+genome
        if tokens[len(tokens)-1] == "gz": cmd = "mkdir "+genome+"; tar xzvf "+found+" -C "+genome
        if cmd != "":          
            result = os.system(cmd)
            print("Result:\t"+str(result))
            if result != 0:
                raise ValueError("External call caused an error, terminating...")


if __name__ == '__main__':
    print("Starting program...")
    # constructing command line parser
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--genome',action='store',type='string',dest='genome',help='Specify the strand column name or number (zero-based) for the left input file')    
    parser.add_option('--genomeDir',action='store',type='string',dest='genomeDir',help='Specify the directory in which the genomes are stored (default=".")',default=".")    
    (options,args) = parser.parse_args()
    
    if not options.genome: 
        print("Mandatory parameters missing. Program will terminate now.")
        print("\nYour parameter settings:")
        print(options)
        raise SystemExit

    downloadGenome(options.genome, options.genomeDir)
    print("Program successfully terminating....")
