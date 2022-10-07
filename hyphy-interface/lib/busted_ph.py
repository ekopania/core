############################################################
# Functions for the BUSTED-PH model from Hyphy.
# 11.2020
############################################################

import sys, os, json, re, lib.hpcore as hpcore, lib.hpseq as hpseq, lib.hptree as tp
from collections import defaultdict

############################################################

#Copies from Gregg Thomas's paml generator scripts for branch-site tests
def assignTargetClade(tree_str, targs):
# Function that reads a tree into a dictionary and determines the target (i.e., foreground) branch from a set of tip labels

    #print("TREE BEFORE PARSING: ", tree_str);
    tinfo, t, root = tp.treeParse(tree_str);
    # Parse tree into dictionary

    #print(tinfo);

    target_found, repl_nodes = False, [];
    for n in tinfo:
        cur_clade = tp.getClade(n, tinfo);
        # Get the clade that defines each node

        if all(c_nodes in targs for c_nodes in cur_clade):
            repl_nodes.append(n);
            target_found = True;
        # If the clade matches the specified targets, set as a target branch

    if target_found:
        for n in repl_nodes:
            t = t.replace(n, n + "{TEST}");
        t = re.sub('<[\d]+>', '', t) + ";";
    # Add target labels and remove node labels

    return target_found, t;

# Same as above, but for reference (i.e., background) branches
def assignReferenceClade(tree_str, refs):

    tinfo, t, root = tp.treeParse(tree_str);
    # Parse tree into dictionary

    ref_found, repl_nodes = False, [];
    for n in tinfo:
        cur_clade = tp.getClade(n, tinfo);
        # Get the clade that defines each node

        if all(c_nodes in refs for c_nodes in cur_clade):
            repl_nodes.append(n);
            ref_found = True;
        # If the clade matches the specified reference, set as a reference branch

    if ref_found:
        for n in repl_nodes:
            t = t.replace(n, n + "{REFERENCE}");
        t = re.sub('<[\d]+>', '', t) + ";";
    # Add reference labels and remove node labels

    return ref_found, t;


############################################################


def generate(indir, tree_input, gt_opt, aln_id_delim, targets, references, hyphy_path, outdir, logdir, outfile):
    if aln_id_delim:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : f.split(aln_id_delim)[0], "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    else:
        aligns = { os.path.splitext(f)[0] : { "aln-file" : os.path.join(indir, f), "id" : "NA", "tree" : False } for f in os.listdir(indir) if f.endswith(".fa") };
    # Read and sort the alignment file names

    for aln in aligns:
        if gt_opt:
            if aln_id_delim:
                tree_dir = os.path.join(tree_input, aln);
                if os.path.isdir(tree_dir):
                    tree_dir_files = os.listdir(tree_dir);
                    tree_file = "".join([ f for f in tree_dir_files if re.findall(aligns[aln]['id'] + '(.*).treefile', f) != [] and "rooted" not in f ]);
                    tree_file = os.path.join(tree_dir, tree_file);
                else:
                    tree_file = False;
            else:
                tree_file = os.path.join(tree_input, aln, aln + ".treefile");
        else:
            tree_file = tree_input;

        if os.path.isfile(tree_file):
            aligns[aln]['tree'] = tree_file;
        # Read the tree and remove any bootstrap node labels.1
    # Read the appropriate tree depending on the -tree and -genetree options.

    tree_skipped, target_skipped, reference_skipped, stop_skipped = 0, 0, 0, 0;
    for aln in aligns:
        if not aligns[aln]['tree']:
            outfile.write(" # Tree file not found. Skipping: " + aln + "\n");
            tree_skipped += 1;
            continue;
        # Check the tree file.          
        
        #print(" # Treefile:" + aligns[aln]['tree']);
        g_tree = open(aligns[aln]['tree'], "r").read().strip();
        g_tree = re.sub("\)[\de.-]+:", "):", g_tree);
        # Read in the tree from treefile.

        target_found, targ_tree = assignTargetClade(g_tree, targets);
        if not target_found:
            outfile.write(" # Target clade not found. Skipping: " + aligns[aln]['aln-file'] + "\n");
            target_skipped += 1;
            continue;
        # Assign the target lineages to the current tree.

        print("#TREE WITH TEST BRANCHES LABELED:")
        print(targ_tree);
        print("#END LABELED TEST TREE:")

        ref_found, cur_tree = assignReferenceClade(targ_tree, references);
        if not ref_found:
            outfile.write(" # Reference clade not found. Skipping: " + aligns[aln]['aln-file'] + "\n");
            reference_skipped += 1;
            continue;
        # Assign the reference lineages to the current tree.

        seq_dict = hpseq.fastaGetDict(aligns[aln]['aln-file']);
        prem_stop_flag = False
        for title in seq_dict:
            prem_stop, new_seq = hpseq.premStopCheck(seq_dict[title], allowlastcodon=True, rmlast=True);
            if prem_stop:
                prem_stop_flag = True;
            seq_dict[title] = new_seq;
        # Read the sequence and check for premature stop codons.

        if prem_stop_flag:
            outfile.write(" # Premature stop found. Skipping: " + aln + "\n");
            stop_skipped += 1;
            continue;
        # Check the alignment for premature stop codons (which PAML hangs on)

        cur_outdir = os.path.join(outdir, aln);
        if not os.path.isdir(cur_outdir):
            os.system("mkdir " + cur_outdir);
        # Make the output directory for this alignment

        new_treefile = os.path.join(cur_outdir, "busted-ph.tre");
        with open(new_treefile, "w") as treefile:
            treefile.write(cur_tree);
        # Make the tree file for this alignment

        new_seqfile = os.path.join(cur_outdir, "busted-ph.fa");
        with open(new_seqfile, "w") as seqfile:
            for title in seq_dict:
                tip_name = str(title).split(">")[1]
                #print(tip_name)
                split_tree = re.split(' |\(|\)|,', cur_tree)
                #print(split_tree)
                if tip_name in split_tree:
                    #print("Writing to output alignment fasta")
                    seqfile.write(title + "\n");
                    seqfile.write(seq_dict[title] + "\n");
        # Write the sequences for this alignment

        cur_jsonfile = os.path.join(outdir, aln + ".json");
        #cur_outfile = os.path.join(cur_outdir, align + "-out.txt");
        cur_logfile = os.path.join(logdir, aln + ".log");
        # Get the output json and log file names

        hyphy_cmd = "hyphy " + hyphy_path + "BUSTED-PH/BUSTED-PH.bf --alignment " + new_seqfile + " --tree " +  new_treefile + " --srv No --branches TEST --comparison REFERENCE --output " + cur_jsonfile + " &> " + cur_logfile 
        outfile.write(hyphy_cmd + "\n");
        # Construct and write the hyphy command

    hpcore.PWS("# Num skipped because tree file not found     : " + str(tree_skipped), outfile);
    hpcore.PWS("# Num skipped because target clade not found  : " + str(target_skipped), outfile);
    hpcore.PWS("# Num skipped because reference clade not found  : " + str(reference_skipped), outfile);
    hpcore.PWS("# Num skipped because of premature stop codons: " + str(stop_skipped), outfile);

############################################################

def parse(indir, features, outfile, pad):

    if features:
        headers = ["file","id","chr","start","end","lrt","pval"];
    else:
        headers = ["file","branch","lrt","pval"];
    outfile.write(",".join(headers) + "\n");
    # Write the output headers 

    hyphy_files = os.listdir(indir);
    num_files = len(hyphy_files);
    num_files_str = str(num_files);
    num_files_len = len(num_files_str);
    # Read align file names from input directory

    num_unfinished = 0;
    # A count of the number of unfinished hyphy runs as determined by empty output files

    counter = 0;
    for f in os.listdir(indir):
        if counter % 500 == 0:
            counter_str = str(counter);
            while len(counter_str) != num_files_len:
                counter_str = "0" + counter_str;
            print ("> " + hpcore.getDateTime() + " " + counter_str + " / " + num_files_str);
        counter += 1;
        # Track progress

        cur_json_file = os.path.join(indir, f);
        if not os.path.isfile(cur_json_file) or not cur_json_file.endswith(".json"):
            continue;
        if os.stat(cur_json_file).st_size == 0:
            num_unfinished +=1 ;
            continue;
        # Get the current output file.

        #print(f);
        if features:
            if "-" in f:
                fid = f.split("-")[0];
            elif "." in f:
                fid = f.split(".")[0];
            else:
                fid = f;
            cur_feature = features[fid];
            # Look up transcript/gene info for current gene to save in output, if metadata is provided.

        with open(cur_json_file) as json_data:
            cur_data = json.load(json_data);
        # Read the Hyphy json file.

        #print(cur_data)

        if features:
            gene_info = { 'id' : fid, 'chr' : cur_feature['chrome'], 'start' : cur_feature['start'], 'end' : cur_feature['end'],
                "lrt" : "NA", "pval" : "NA" };
        else:
            gene_info = { "lrt" : "NA", "pval" : "NA" };   
        # Initialize the output dictionary for the current branch.

        #gene_info["dn/ds"] = str(cur_data["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]);
        gene_info["lrt"] = str(cur_data["test results"]["LRT"]);
        gene_info["pval"] = str(cur_data["test results"]["p-value"]);
        # Retrieve the rate estimates from the json data.

        gene_outline = [f] + [ gene_info[h] for h in headers if h not in ["file"] ];
        outfile.write(",".join(gene_outline) + "\n");
        # Ouput rate estimates for both the gene.

    hpcore.PWS("# ----------------", outfile);
    #hpcore.PWS(hpcore.spacedOut("# Number unfinished:", pad) + str(num_unfinished), outfile);
