# DUSmapper

DUSmapper is a Python script which identifies putative genomic islands in genome sequence data (primarily for Neisseriaceae species) via atypical DNA uptake sequence distributions.
It is described further in [this paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000372) (please cite this paper if you use DUSmapper in your publications).   

Depending on your preferred way of using python I have provided a .py script to be run from the command line, and an alternative .ipynb jupyter notebook.  

Tested on a Mac with python 3.9. 


## DEPENDENCIES

The following prerequisite python packages are included if you install python via [anaconda](https://www.anaconda.com/):

- Biopython (https://biopython.org/)
- Numpy
- Matplotlib

Whereas the following need to be installed separately:  

- fuzznuc (from the EMBOSS package) – a pattern search tool for nucleotide sequences (http://emboss.toulouse.inra.fr/cgi-bin/emboss/help/fuzznuc)
- DNAPlotter – used to visualise circular genome plots of the script output (https://www.sanger.ac.uk/science/tools/dnaplotter)

Before running the script it is a good idea to make sure fuzznuc is running okay in the command line (I installed fuzznuc via anaconda 
using `conda install -c bioconda emboss`). **Please note, fuzznuc will not tolerate blank spaces in the path or file names of your genome sequence files, so be sure to replace these with e.g. underscores.** Also try opening some genome sequence files in DNAPlotter and experiment with the track manager.


## SCRIPT INFORMATION

Thes script takes a directory containing single-entry Genbank sequence files as input.
All ocurrences of the specified DUS dialect sequences (allowances for nucleotide mismatches are optional) are identified and
added into the Genbank file as annotations. By default, the specified dialect sequences are the eight 12 bp consensus sequences outlined 
in [this paper](https://doi.org/10.1371/journal.pgen.1003458), but these can be modified.

The most prevalent DUS dialect type (a.k.a primary DUS or prDUS, which is presumably the native, endogenous
DUS targeted by ComP of the species) is identified, and regions of the genome containing a low density of prDUS ('DUS islands')
are also annotated. The currently implemented method for desgnating DUS islands is very basic – only the gaps between directly adjacent 
DUS on each strand are calculated, rather than the true local DUS density for a windowed region. DUS gaps which are greater than the
median gap spacing + 4*IQR (interquartile range) are classified as DUS islands, but this cutoff was chosen empirically by comparison with previously characterised 
Neisseria genomic islands, rather than by rigorous statistics. 
I plan to implement a better method for identifying DUS islands based on local frequencies in future.

In addition to annotating DUS dialects and prDUS islands, the script enables visualision of the resulting annotated genome
with DNAPlotter by outputting a template file (specifying the number of tracks and their colours etc.) and a shell script 
for opening the annotated genome file in DNAPlotter. DNAPlotter is Java-based, and your Java location may need to be added 
to the system path. In order to generate the DUSPlotter template file and shell launcher script, the path to the blank template 
file and path to your DNAPlotter application .jar file must be provided (either as the 2nd and 3rd command line arguments, respectively, for the
.py file, or by modifying the path strings in the .ipynb notebook). A blank template file 'blank8DUStemplate' is provided, but you will need
to generate your own custom template if you wish to plot different features, or different DUS dialects etc.
You may prefer to omit this DNA plotter output and instead visualise genomic DUSs and DUS islands using different software.

Several plots of prDUS gap distributions are also saved, as well as .csv files listing counts for each DUS dialect and some prDUS
gap statistics.

Despite a few bells and whistles, at its heart this script just searches for DUS motifs within a genome sequence file
(using fuzznuc), detects prDUS islands, and annotates the sequence file with these features.


## RUNNING THE CODE

DNAPlotter.py can be run from the command line with only a single argument – the directory containing each genome to be analysed as a single-entry
genbank file (but in this case no DUSPlotter template or launcher will be output), e.g.:

`DUSplotter.py /path/to/genome/directory`

To output the template file and launcher bash script for DNAPlotter, two additional arguments must be provided after the genome directory path: 
(1) the path to the blank DNAPlotter template file, and (2) the path to your DNAPlotter application .jar file, e.g.:

`DUSplotter.py /path/to/genome/directory /path/to/blank8DUStemplate /path/to/DNAPlotter.app/Contents/Java/DNAPlotter.jar`

If you opt to run the jupyter notebook version of the script, you just need to assign the relevant paths to the `dir_path`, `plotter_template`, 
and `DNAPlotter_path` variables.

I have provided three example genomes (Neisseria meningitidis MC58, Neisseria lactamica 020-06, and Wielerella bovis CCUG 44465) to test the script with, as well as a folder containing the expected script output from each genome. In the 'expected output' folders, the launch scripts and templates for DUSPlotter won't actually function because they contain absolute paths to the specific genome files I used when running DUSPlotter.

## NOTES

If you have multiple-entry GenBank genome files, these need to be converted to single-entry files before running DUSplotter. I have provided a short 
script ('Multi -> single-record Genbank converter') for doing this. The script takes the directory containing your genome files as input, and outputs a 
new directory containing concatenated single-entry versions (by default, entries are concatenated with 50 intervening 'N' characters to clearly denote
contig boundaries).

Assignment of the 'plus' and 'minus' strand by the script is arbitrary, and is not based on e.g. replication or transcription directions.

The `blank8DUStemplate` DUSPlotter template file specifies the following tracks from innermost to outermost: 

(1) GC skew  
(2) GC content  
(3) TG-wadDUS – Grey  
(4) AA-king3DUS – Yellow/Grey  
(5) AG-eikDUS – Fuchsia  
(6) AG-mucDUS – Fluorescent Green  
(7) AG-simDUS – Purple  
(8) AG-kingDUS – Orange  
(9) AG-DUS – Blue  
(10) AT-DUS – Green  
(11) CDS forward – Red  
(12) CDS reverse - Red  
(13) DUS islands plus strand – Turquoise  
(14) DUS islands minus strand – Grey/Turquoise  

DNAPlotter is quite sensitive to formatting within Genbank files, and some annotations may need to be edited before the file will open okay in DNAPlotter. For example, forward slashes within feature qualifiers tend to be problematic and should be replaced with e.g. underscores.

If you have difficulties opening the annotated Genbank file in DNAPlotter, just use any other genome plotting software compatible with Genbank files.
