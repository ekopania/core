#!/usr/bin/python
########################################################################################
# A script to do many Newick tree related things
#
# Dependencies: core
#
# Gregg Thomas, Summer 2017
########################################################################################

import sys, os, random, argparse
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/corelib/"))
import core, treeparse as tp, treelib as tree

####################
parser = argparse.ArgumentParser(description="A script to do many Newick tree related things");
parser.add_argument("-i", dest="input", help="A directory containing many Newick tree files or a single file containing Newick tree(s).", default=False);
parser.add_argument("-o", dest="output", help="Desired output location. If input is a file this should be a file, if a directory this will be a directory.", default=False);

parser.add_argument("--sep", dest="tree_sep", help="Given a single file with many Newick trees (one per line), this will write them all to separate files in the output directory.", action="store_true");
parser.add_argument("--join", dest="tree_join", help="Given an input directory with many Newick tree files, this will combine them all into a single file.", action="store_true");
parser.add_argument("--label", dest="label_tree", help="Given an input file or tree string, this will add internal node labels to the tree(s).", action="store_true");
parser.add_argument("--rootcheck", dest="root_check", help="Given an input file or tree string, this will check if the tree(s) are rooted or not.", action="store_true");
parser.add_argument("--root", dest="root_tree", help="Given an input file or tree string, this will root the tree with the specified outgroup(s).", action="store_true");
parser.add_argument("--concordance", dest="fotc", help="Given an input species tree and a file containing many single-copy gene trees this module will calculate concordance factors for each node in the species tree. Use -genetrees for the input gene tree file and -i for the input species tree file or string.", action="store_true");
parser.add_argument("--tipcount", dest="count_tips", help="Given a file with many trees, simply count the number of unique tip labels in all trees.", action="store_true");
parser.add_argument("--relabeltips", dest="relabel", help="Given a file with many trees and a set of labels defined by -labels, this will relabel tip nodes.", action="store_true");

parser.add_argument("-prefix", dest="file_prefix", help="For --sep, a string that will be used as the base file name for each output file.", default=False);
parser.add_argument("-outgroup", dest="outgroup", help="For --root, a comma separated list of tip labels common between trees to use as the outgroup for rooting", default=False);
parser.add_argument("-genetrees", dest="genetrees", help="For --concordance, this is the file containing the gene trees.", default=False);
parser.add_argument("-labels", dest="labels", help="For --relabeltip, the old label and the newlabel in the format: \"old1,new1 old2,new2\". Old labels don't need to match exactly with existing labels to allow for matching substrings.", default=False);

args = parser.parse_args();
# Input option definitions.

if not args.input or not os.path.exists(args.input):
	if args.label_tree or args.root_check or args.root_tree or args.fotc:
		file_flag = False;
		tree_flag = True;
	elif not os.path.exists(args.input):
		sys.exit(core.errorOut(1, "-i must be specified and must be a valid file or directory name."));	
else:
	if os.path.isfile(args.input):
		file_flag = True;
		tree_flag = False;
		filelist = [os.path.abspath(args.input)];
	else:
		file_flag = False;
		tree_flag = False;
		filelist_init = os.listdir(args.input);
		filelist = [os.path.abspath(os.path.join(args.input, f)) for f in filelist_init];
# This checks if the input (-i) entered is valid. If so, it parses it as either a directory or a single file.

if not args.output and not tree_flag:
	print "\n** Warning -- No output location specified. Will be determined automatically.";

pad = 40;
if args.tree_sep:
	if not file_flag:
		sys.exit(core.errorOut(2, "--sep only works on an input FILE!"));
	output, outnum = core.defaultOutDir(args.input, file_flag, "sep", args.output);
	print "* Making output directory...";
	os.system("mkdir " + output);
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Separating trees in:", pad), args.input;
	print core.spacedOut("Writing output to:", pad), output;
	print "-------------------------";
	tree.treeSep(filelist[0], args.file_prefix, output);
	sys.exit();
# --sep : takes an input FILE with many trees and puts each tree in its own file.

if args.tree_join:
	if file_flag:
		sys.exit(core.errorOut(3, "--join only works on an input DIRECTORY containing many tree files!"));
	output, outnum = core.defaultOutFile(args.input, file_flag, "join", args.output);
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Joining trees in:", pad), args.input;
	print core.spacedOut("Writing output to:", pad), output;
	print "-------------------------";
	tree.treeJoin(filelist, output);
	sys.exit();
# --join : takes an input DIRECTORY with many files and combines all trees from all files into a single file.

if args.label_tree:
	if file_flag == False and tree_flag == False:
		sys.exit(core.errorOut(4, "--label only works on an input FILE containing many trees or a TREE STRING."));
	if tree_flag:
		try:
			td, t, root = tp.treeParse(args.input);
			filelist = [td,t,root];
			# Weird hack.
			output = "";
		except:
			sys.exit(core.errorOut(5, "Couldn't read the input string as a Newick tree!"));
	# If the input is not a file, check to see if it's a Newick string.
	else:
		output, outnum = core.defaultOutFile(args.input, file_flag, "label", args.output);
		print "=======================================================================";
		print "\t\t\t" + core.getDateTime();
		print core.spacedOut("Labeling all trees in:", pad), args.input;
		print core.spacedOut("Writing labeled trees to:", pad), output;
	tree.labelTree(filelist, tree_flag, output);
	sys.exit();
# --label : takes an input Newick string or file and puts labels on the internal nodes.

if args.root_check:
	print "\n** Warning -- The root check module works on the basis that an unrooted tree has a trifurcation at the 'root.' So if your tree has other non-bifurcating nodes this will not give reliable results!";
	print "\n              Only use this module with bifurcating trees!";
	if file_flag == False and tree_flag == False:
		sys.exit(core.errorOut(6, "--rootcheck only works on an input FILE containing many trees or a TREE STRING."));
	if tree_flag:
		try:
			td, t, root = tp.treeParse(args.input);
			filelist = [td,t,root];
			# Weird hack.
			output = "";
		except:
			sys.exit(core.errorOut(7, "Couldn't read the input string as a Newick tree!"));
	# If the input is not a file, check to see if it's a Newick string.
	else:
		output, outnum = core.defaultOutFile(args.input, file_flag, "rootcheck", args.output);
		print "=======================================================================";
		print "\t\t\t" + core.getDateTime();
		print core.spacedOut("Checking all trees in:", pad), args.input;
		print core.spacedOut("Writing output to:", pad), output;
	tree.rootCheck(filelist, tree_flag, output);
	sys.exit();
# --rootcheck : takes an input Newick string or file and checks if the trees are rooted or not.

if args.root_tree:
	print "\n** Warning -- The root option relies on the software Newick Utilities. If you don't have this installed and in your PATH variable, you will see an error!";
	if file_flag == False and tree_flag == False:
		sys.exit(core.errorOut(8, "--root only works on an input FILE containing many trees or a TREE STRING."));
	if not args.outgroup:
		sys.exit(core.errorOut(9, "With --root, a set of tip labels must be specified as the desired outgroup (-outgroup)."));
	if tree_flag:
		try:
			td, t, root = tp.treeParse(args.input);
			filelist = [td,t,root,args.input];
			# Weird hack.
			output = "";
		except:
			sys.exit(core.errorOut(10, "Couldn't read the input string as a Newick tree!"));
	# If the input is not a file, check to see if it's a Newick string.
	else:
		output, outnum = core.defaultOutFile(args.input, file_flag, "reroot", args.output);
		print "=======================================================================";
		print "\t\t\t" + core.getDateTime();
		print core.spacedOut("Re-rooting all trees in:", pad), args.input;
		print core.spacedOut("Writing rooted trees to:", pad), output;
	tree.rootTrees(filelist, tree_flag, args.outgroup, output);
	sys.exit();
# --root : takes an input Newick string or file and roots or re-roots the trees using Newick Utilities and the specified outgroups.

if args.fotc:
	if file_flag == False and tree_flag == False:
		sys.exit(core.errorOut(11, "--concordance only works on an input FILE containing many trees or a TREE STRING."));
	if not os.path.isfile(args.genetrees):
		sys.exit(core.errorOut(12, "-genetrees must be a valide file name!"));
	else:
		args.genetrees = os.path.abspath(args.genetrees);
	# Check if the input files are valid.

	if tree_flag:
		try:
			td, t, root = tp.treeParse(args.input);
			filelist = [td,t,root,args.input];
			# Weird hack.
		except:
			sys.exit(core.errorOut(13, "Couldn't read the input string as a Newick tree!"));
	# If the input is not a file, check to see if it's a Newick string.
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print "Calculating concordance factors for your species tree.";
	print core.spacedOut("Using gene trees in:", pad), args.genetrees;
	print "Simply printing output to the screen";
	tree.flightOfTheConcordance(filelist, tree_flag, args.genetrees);
	sys.exit();
# --concordance : takes an input species tree (Newick string or file) and single-copy gene trees (file) 
# and calculates concordance factors for each internal node of the species tree.

if args.count_tips:
	if not file_flag:
		sys.exit(core.errorOut(14, "--tipcount takes an input (-i) FILE only."));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Counting tips in:", pad), args.input;
	tree.countTips(filelist[0]);
	sys.exit();
# --tipcount : counts all unique tips in an input file with many trees.

if args.relabel:
	if not file_flag:
		sys.exit(core.errorOut(14, "--relabeltips takes an input (-i) FILE only."));
	if not args.labels:
		sys.exit(core.errorOut(15, "-labels must be entered with --relabeltips"));
	print "=======================================================================";
	print "\t\t\t" + core.getDateTime();
	print core.spacedOut("Relabeling tips in:", pad), args.input;
	output, outnum = core.defaultOutFile(args.input, file_flag, "relabel", args.output);
	print core.spacedOut("Writing output to:", pad), output;
	tree.relabelTips(filelist[0], args.labels, output);
	sys.exit();
# --relabeltips : in a file containing many trees, relabels all tips containing an old label with a new label specified by user
