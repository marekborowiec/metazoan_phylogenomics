#! /usr/bin/env python

# phylo_jackknife.py by Marek Borowiec
# 15 December 2014

# This program draws and concatenates random samples 
# from a pool of single gene alignments
# for a phylogenetic jackknife analysis

import glob, shutil, random, os, subprocess

try:
    import readline
except ImportError:
    print "Module readline not available."
else:
    import rlcompleter

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

readline.set_completer_delims('\t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

prompt = "> "

print("""What is the path to your working directory?
It should contain only single-locus alignments in FASTA or NEXUS format. 
example: /home/user/Data/AlignedFiles""")
dirPath = raw_input(prompt)

print("""How many loci do you want to use for your replicate?
example: 20""")
lociNo = raw_input(prompt)

print("""How many jackknife replicates do you want to perform?
example: 200""")
replicatesNo = raw_input(prompt)

print("""What is the path to your phyutility.jar file?
example: /home/user/Phylo-Software/Phyutility/phyutility.jar""")
phyutDir = "java -jar " + raw_input(prompt) + " -concat -in "
#print phyutDir

print("""Are your sequences amino acid or DNA?
(aa/dna)""")
seqType = raw_input(prompt)


filenames = [f for f in os.listdir(dirPath) if os.path.isfile(dirPath + '/' + f)]
#print(filenames)
#print(len(filenames))
for  i in range(1, int(replicatesNo) + 1):

	randomFiles = random.sample(filenames, int(lociNo))
	#print(randomFiles)
	destDirectory = dirPath + '/Sample_' + lociNo + '_' + str(i) + '/'
	print(destDirectory)
	os.mkdir(destDirectory)

	for fname in randomFiles:

    		srcPath = os.path.join(dirPath, fname)
		destFile = destDirectory + fname
    		shutil.copyfile(srcPath, destFile)

directories = [x[0] for x in os.walk(dirPath)]
print directories

for directory in directories:
        #print directory
        call_string = phyutDir + directory + '/* ' + '-out ' + directory + '.nex'
        print call_string
        subprocess.call(call_string, shell=True)
