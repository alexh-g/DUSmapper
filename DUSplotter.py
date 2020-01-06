from Bio.Emboss.Applications import FuzznucCommandline
from Bio import SeqIO
import os, sys
from os.path import basename, splitext, join
import glob
from shutil import copyfile

#search for DUSs and save report as genbank feature file
DUSs = [['AT-DUS','ATGCCGTCTGAA'],["AG-DUS", "AGGCCGTCTGAA"],["AG-kingDUS", "AGGCAGCCTGAA"],
        ["AG-mucDUS", "AGGTCGTCTGAA"],["AG-simDUS", "AGGCTGCCTGAA"],
        ["AG-eikDUS", "AGGCTACCTGAA"],["AA-king3DUS", "AAGCAGCCTGCA"],
        ["TG-wadDUS", "TGCCTGTCTGAA"]]

for filename in glob.glob(join(sys.argv[1],'*.gb*')):
    genome_base = splitext(basename(filename))[0]
    outdir = join(os.path.dirname(filename),genome_base+"_fuzz_output")
    if not os.path.exists(outdir):
        os.mkdir(outdir)

        for i, DUS in enumerate(DUSs):
            fuzz = FuzznucCommandline()
            fuzz.sequence = filename
            fuzz.complement = True
            fuzz.pattern = DUS[1]
            fuzz.pmismatch = 0
            fuzz.rformat = "genbank"
            outfile = join(outdir, DUS[0]+"_fuzz_out.gbk")
            fuzz.outfile = outfile
            std_out, err_output = fuzz()

            #edit notes of saved file:
            edit = 'fuzz_out_'+DUS[0]+'.gbk'
            with open(outfile, 'r') as fuzz_results, open(join(outdir, edit), 'w') as w:
                for line in fuzz_results:
                    if not 'nucleotide_motif' in line:
                        line = line.replace('*pat pattern:'+DUS[1], DUS[0])
                        w.write(line)
                    else:
                        #add gene field so DUS type is displayed in Artemis
                        line = line.replace("/note=\"*Type SO:0000714 nucleotide_motif\"", "/gene=\""+DUS[0]+"\"")
                        #print(line) #check output
                        w.write(line)
            os.remove(outfile) #delete clunkily annotated output files

            #paste DUSs into genome:
            if DUS[0]=='AT-DUS':
                for_annotation = join(outdir, genome_base+"_for_"+DUS[0]+"_annotation.gb")
                copyfile(filename, for_annotation)
                pastefile = for_annotation
            else:
                pastefile = join(outdir,genome_base+DUSs[i-1][0]+'_annotated.gb')

            if DUS[0]=='TG-wadDUS':
                out = join(outdir,genome_base+'_annotated.gb')
            else:
                out = join(outdir,genome_base+DUS[0]+'_annotated.gb')

            with open(pastefile, 'r') as genome, open(out, 'w') as w:
                for line in genome:
                    #remove any plasmids from multi-record genbank files by breaking at the // line
                    #in rare cases the plasmid sequence precedes the gDNA, so this method doesn't work
                    if "//\n" not in line:
                        if not 'Location/Qualifiers' in line:
                            w.write(line)
                        else:
                            DUS_file = open(join(outdir, edit), 'r')
                            #remove first FEATURES line from fuzznuc output
                            DUS_feat = ''.join(l for l in DUS_file.readlines() if 'Qualifiers' not in l)
                            line += DUS_feat
                            w.write(line)
                    else:
                        w.write(line)
                        break

            os.remove(pastefile)
            os.remove(join(outdir, edit))

        #create bespoke template file :
        original_template = '/path/to/blank8DUStemplate' #replace with file path DUSPlotter template file
        template_edit = join(outdir, genome_base+'8DUStemplate')
        seq = SeqIO.read(out,'genbank')

        with open (original_template, 'r') as r, open(template_edit,'w') as w:
            for line in r:
                if 'end=' in line:
                    line = line.replace('length', str(len(seq)))
                    w.write(line)
                elif not 'genome_file' in line:
                    w.write(line)
                else:
                    line = line.replace('genome_file', genome_base+'_annotated.gb')
                    line = line.replace('path', os.path.abspath(outdir))
                    w.write(line)

        #create bespoke DNAPlotter launcher:
        launcher = '/path/to/DNAplotter_launch_template.sh'
        launch_out = join(outdir, genome_base+'plot.command')

        with open(launcher, 'r') as r, open(launch_out,'w') as w:
            for line in r:
                line = line.replace('insert_path', os.path.abspath(template_edit))
                w.write(line)
        #make bash script executable:
        os.system('chmod 755 '+os.path.abspath(launch_out))
