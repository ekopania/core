# CORE: Code fOr Romps in Evolutionary data
### A mixture of scripts and libraries to help with sequence data manipulation, tree parsing, and other things.

## Author
#### Gregg Thomas

### About
###### These scripts can be used for many tasks including sequence handling, tree making, and sequence alignment.
###### Many of these programs are mainly used to easily run other genomics or phylogenetics programs on a bunch of files. Pay attention to the dependencies for each script to make sure you have the proper programs installed.
###### Please note that many of these scripts expect input as FASTA files. For my scripts, these *must* have the extension .fa. If you don't have FASTA formatted files, you can use seq_convert to get them to FASTA format and fa_edit to make any changes you need to them afterwards.

###### Almost all of these scripts are written in Python 2.7 (https://www.python.org/downloads/).
###### For any script, use the -h flag for specific usage details.

1. corelib/core.py
  * General helper functions such as reading sequences to a dictionary. You'll have to look to see what all is there.
2. corelib/nj_tree.r
  * Simple R script to get a Neighbor Joining tree. Used by supertreemaker and probably not helpful standalone.
  * Dependencies:
  <ol>
   <li>R (https://www.r-project.org/)</li>
  </ol>
3. corelib/treeparse.py
  * A couple functions that read newick formatted trees and return all relevant information in a more useful way to code with.
4. count_aln.py
  * This script gathers statistics about a single alignment file, or a directory full of alignment files.
5. count_pos.py
  * This script simply counts the number of amino acids or nucleotides in a file or directory.
6. fa_concat.py
  * Concatenates many FASTA formatted sequence files into a single FASTA file.
7. fa_edit.py
  * A general purpose FASTA handling script. Can relabel and trim headers and remove start and stop AAs.
8. get_bs.py
  * A script to get the bootstrap replicates from a run_raxml run.
9. how\_many\_trees
  * Just a little script to show the number of possible rooted tree topologies for a given number of species.
10. run_codeml.py
  * A script to run some basic PAML analyses with codeml. Still adding a lot to this one. You'll need codeml and newickutils in your PATH.
  * Dependencies:
  <ol>
   <li>PAML (http://abacus.gene.ucl.ac.uk/software/paml.html) called as `codeml`</li>
   <li>Newick Utilities (http://cegg.unige.ch/newick_utils) called as `nw_prune`</li>
  </ol>
11. run_gblocks.py
  * A script to run GBlocks to mask a directory full of alignments in FASTA format. Note: This currently runs GBlocks at the most relaxed settings for phylogenetic tree inference. It will reject any masks that remove more than 20% of the columns from the original alignment.
  * Dependencies:
  <ol> 
   <li>GBlocks (http://molevol.cmima.csic.es/castresana/Gblocks.html) called as `gblocks`</li>
  </ol>
12. run_muscle.py
  * This will make MUSCLE alignments out of a directory of FASTA files.
  * Dependencies:
  <ol>
   <li>muscle (http://www.drive5.com/muscle/downloads.htm) called as `muscle`</li>
  </ol>
13. run\_pasta\_aln.py	
  * This will make PASTA alignments out of a directory of FASTA files.
  * Dependencies:
  <ol>
   <li>PASTA (http://www.cs.utexas.edu/~phylo/software/pasta/) called as `python run_pasta.py`</li>
  </ol>
14. run_raxml.py
  * Runs some basic RAxML analyses on a directory full of FASTA files.
  * Dependenceies:
  <ol>
   <li>RAxML (http://sco.h-its.org/exelixis/web/software/raxml/index.html) called as `raxml`, though you can specify the path to your own raxml executable.</li>
  </ol>
15. seq_convert.py
  * A sequence file format conversion tool. Currently converts between FASTA (.fa), Phylip (.ph), and Nexus (.nex) formats. It assumes files will have those extensions. Remember, these formats vary a lot in the details, so they might not work right away for everything. Let me know if you run into problems and I'll try to fix it.
16. supertreemaker.py
  * This script can do several things. Runs SDM to get average consensus distance matrices, makes NJ trees from distance matrices, re-roots trees, and makes ultrametric trees. All of these programs will need to be in your PATH.
  * Dependencies:
  <ol>
   <li>SDM (http://www.atgc-montpellier.fr/sdm/) to calculate the matrix, called as `java -jar ~/bin/SDM/SDM.jar`</li>
   <li>R (https://www.r-project.org/) to make NJ trees, called as `Rscript`</li>
   <li>Newick Utilities (http://cegg.unige.ch/newick_utils) to re-root trees, called as `nw_reroot`</li>
   <li>r8s (http://loco.biosci.arizona.edu/r8s/) to smooth trees, called as `r8s`</li>
  </ol>
