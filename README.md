# DUSmapper
A Unix-based Python 2.7 script which identifies putative genomic islands in Neisseriaceae genome sequence data via atypical DNA uptake sequence distributions.

Dependencies:

- Biopython (https://biopython.org/)
- DNAPlotter is required for visualisation of the script output (https://www.sanger.ac.uk/science/tools/dnaplotter)

The script takes a directory containing all the genome sequences to be analysed (in GenBank format with each genome concatenated to a single entry) as its only argument. 

Three output files are produced for each input genome: 
- A copy of the input genome annotated with each occurence of the eight 12 bp Neisseriaceae DUS sequences outlined by Frye et al. [1] (by default only exact matches of DUS sequences are identified, but mismatch allowances can be included by adjusting the relevant fuzznuc argument).
- A DUSPlotter template file for plotting each DUS sequence on a different track.
- A shell script for launching DNAPlotter with the corresponding annotated GenBank sequence and template file (the script includes a chmod command which should ensure this file is executable).


[1] Frye SA, Nilsen M, Tønjum T, Ambur OH. Dialects of the DNA uptake sequence in Neisseriaceae. PLoS Genet 2013;9:e1003458.
