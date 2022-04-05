# SingleNucleotidem6Aseq
This repository contains R scripts used for analysis of single-nucleotide-resolution m6A sequence data, including preprocessing, calculating transcript entropy, extracting annotation information from GTF, running statistical tests, determining gene ontologies and generating unique graphics. 
It was created to partially fulfill the requirements to obtain the degree Bachelor of Science with Honors in Computational Biology from the Center for Computational Molecular Biology at Brown University. 

IMPORTANT NOTE: The exact GTF annotation file used for these analyses must match the one used by the sequencing company that provides you with index location of the m6A sites. Otherwise you will have incomplete subsets of genes you can analyze visually using the makeTxPlots.R and TxPlotFns.R scripts, and all the diagrams that you make will have inaccurate CDS start and end sites, and may create unplottable negative indices. 
