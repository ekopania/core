#!/usr/bin/env python3
########################################################################################
# This script converts GFF or GTF files to a more sensible tab delimited file.
#
# Dependencies: core
#
# Gregg Thomas, Fall 2020
########################################################################################

import sys, os, re, argparse
from collections import defaultdict
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core

########################################################################################
# Functions


########################################################################################
# Main

parser = argparse.ArgumentParser(description="GTF/GFF to generic tab conversion.");
parser.add_argument("-i", dest="input", help="The input GTF or GFF formatted file.", default=False);
parser.add_argument("-o", dest="output", help="The output tab delmited file.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="Set this to indicate you wish to overwrite files specified by -outcsv and -outtxt if they already exist. WARNING: This means the original contents of the file will be deleted.", default=False, action="store_true");
#parser.add_argument("--header", dest="header", help="Set to extract the header info for each file.", default=False, action="store_true");
args = parser.parse_args();
# Input option definitions.

if not args.input or not os.path.isfile(args.input):
    sys.exit(core.errorOut(1, "Cannot find input file (-i)."));

if not args.output:
    sys.exit(core.errorOut(4, "Please specify the name of an output file (-o)."));
elif os.path.isfile(args.output) and not args.overwrite:
   sys.exit(core.errorOut(5, "Output file (-o) already exists! Explicity specify --overwrite to overwrite it."));
# I/O parsing and error checking.

pad = 25;
with open(args.output, "w") as outfile:
    core.runTime("# GFF/GTF feature counting.", outfile);
    core.PWS(core.spacedOut("# Input file:", pad) + args.input, outfile);
    core.PWS(core.spacedOut("# Output file:", pad) + args.output, outfile);
    core.PWS("# ----------------", outfile);

    core.PWS("# " + core.getDateTime() + " Detecting compression...", outfile);
    reader = core.getFileReader(args.input);
    if reader == open:
        line_reader = core.readLine;
        read_mode = "r";
    if reader != open:
        line_reader = core.readGzipLine
        read_mode = "rb";

    core.PWS("# " + core.getDateTime() + " Counting features...", outfile);
    features, total_feature= {}, 0;
    for line in reader(args.input, read_mode):
        #print(line);
        line = line_reader(line);
        if "##FASTA" in line[0]:
            break;
        if line[0][0] == "#":
            continue;
        feature_type, start, end = line[2], int(line[3]), int(line[4]);
        cur_len = end - start;

        if feature_type not in features:
            features[feature_type] = { 'count' : 0, 'lens' : [] };

        features[feature_type]['count'] += 1;
        features[feature_type]['lens'].append(cur_len);
        total_feature += 1;
    core.PWS(core.spacedOut("# Features types read:", pad) + str(len(features)), outfile);
    core.PWS(core.spacedOut("# Total features read:", pad) + str(total_feature), outfile);
    core.PWS("# ----------------", outfile);
    # Read genes.

    core.PWS("# " + core.getDateTime() + " Writing output...", outfile);
    core.PWS("# ----------------", outfile);
    headers = "feature type,count,avg length";
    outfile.write(headers + "\n");
    # Output headers.

    for f in features:
        outline = f + "," + str(features[f]['count']) + "," + str( float(sum(features[f]['lens'])) / float(len(features[f]['lens'])) );
        outfile.write(outline + "\n");
    print("Done!");