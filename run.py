import chimeramate_cy
import sys
args = sys.argv
'str otu_file, int k, int cutoff, k_chimsearch, abskew'
softmasked_file =  args[1].split('.fasta')[0] + '_softmasked.fasta'
d = chimeramate_cy.kmer_filter(args[1],int(args[2]), int(args[3]))
d.softmask(d.reads, d.kmer_abundance_sorting, softmasked_file)

e = chimeramate_cy.chimera_search(softmasked_file, int(args[4]), float(args[5]))