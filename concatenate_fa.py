#!/usr/bin/python
#############################################################################
#This script concatenates FASTA formatted sequences in multiple files to a single
#FASTA file.
#
#Usage: python count_aln.py [input directory] [output file name]
#
#Depenencies: core
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

ins = sys.argv[1];
outs = sys.argv[2];

print "=======================================================================";
print "Concatenating alignments in:\t\t" + ins;
print "Writing concatenated alignment to:\t" + outs;
print "-------------------------------------";

filelist = os.listdir(ins);

concats = {};

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	numbars, donepercent = core.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;

	if each.find(".fa") == -1:
		continue; 	

	infilename = ins + each;

	inseqs = core.fastaGetDict(infilename);

	for title in inseqs:
		newtitle = title[:title.index(" ")];
		if newtitle not in concats:
			concats[newtitle] = inseqs[title];
		else:
			concats[newtitle] = concats[newtitle] + inseqs[title];

outfile = open(outs, "w");
for spec in concats:
	outfile.write(spec);
	outfile.write("\n");
	outfile.write(concats[spec]);
	outfile.write("\n");
outfile.close();

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print "\nDone!";
print "=======================================================================";
