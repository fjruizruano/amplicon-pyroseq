amplicon-pyroseq
================

Processing and analyzing pyrosequencing reads from amplicons

# Dependencies

* EMBOSS
* Acacia
* Uchime
* MAFFT

# Pipeline


We wrote a series of custom Python scripts and used them for counting the number of reads corresponding to each sequence (haplotype) in each sample. This was done in five consecutive steps that we recommend to run inside different folders:

1) To select the ITS2 amplicons, local alignments of the initial half of the sequencing reads against the forward and reverse primers were carried out using the Smith-Waterman algorithm, as implemented in the EMBOSS suite (Rice et al., 2000). We only considered reads with high identity with one of the primers. Additionally, this alignment served to determine the orientation of each read in respect of the aligned primer, since we had reads in both orientations because the 454 adaptors were not included in the primers. We substituted the reads matched to the reverse primer with their reverse-complementary in order to have all the amplicons in the same orientation. Before running, add the sequence of your forward and reverse primers (for ITS2, primers its3 and its4) in lines 65-72 of and the name of your FASTQ file in lines 80 and 97 of the first script.
```
$ python 01split_genes.py
```
The output will be several FASTQ files. "noMat.fq" corresponds to the reads without a good match with the introduced primers. In the other ones you will found the sequences matched with each pair of primers.
2) We assigned each read to a given sample according to its tag combination using local alignments of the complete reads againts the tagged primers. 6mer Tags were designed using using EDITTAG (Faircloth and Glenn, 2012), differing in four or more nucleotides. To avoid missassignation, we sorted the alignments from higher to lower similarity and assigned the sample to a given read only if the best alignment showed less than four sequence differences, (i.e. the minimal number of differences between tags) in both tagged primers. Selected reads were stored in different files according to their combination of tags. The alignments also allowed us to trim tags and primers from the sequences. Before running the script, check the name of the single region you want to analyze and add the tagged primers used in the experiment (for ITS2, 79-97 lines). 
```
$ python 02split_tag.py
```
The program will ask the name of the gene ("its2" for the ITS2 region) and the allowed number of differences to accept a match. The output is a FASTQ file for each tag combination, i. e. sample.
3) We used the Acacia software (Bragg et al., 2012) to correct typical 454 pyrosequencing errors, i.e. substitutions and indels in homopolymeric regions. We then generated a file with all found haplotypes per sample and their frequency, and we searched for chimeric sequences with UCHIME (Edgar et al., 2011) using the default options for the de novo algorithm. To run this script, create a "lista.txt" file in the working directory with the path of each generated file in the previous step. It can be created with "$ ls out_its3* > lista.txt".
```
$ python 03clean_amp.py
```
It will generate for each sample a "*_hap.fas" FASTA file with all the found haplotypes and their abundance in the UCHIME format and a "*_hap_chime.out" file as output of UCHIME to be manually revised.
4) In the chase of the ITS2, toannotate the 5.8S, ITS2 and 28S regions in all the haplotype files, per sample, we aligned the sequences using MAFFT v7 (Katoh and Standley, 2013) with LINSI options including an additional E. plorans ITS2 sequence (GenBank accession number JN811835.1, ref.fas file) as a reference, annotating its ITS2 region with the ITS2 Database III annotation tool (Koetschan et al., 2010). Since the 454 reads included partial regions of the 5.8S and 28S rRNA genes, in addition to the ITS2 sequence, we used these partial coding regions as internal control for sequencing errors. This provided a way to avoid false positives in identifying genuine ITS2 haplotypes (i.e. variants for the ITS2 region sequence) in the gDNA. For stringency, we considered as sequencing errors all the variation found in these region, compared to the sequence JN811835.1. The proportion of reads carrying any variation in respect to the conserved coding sequence was thus our estimate of the maximum error rate of the experiment. In order to select the genuine ITS2 haplotypes in the different samples of the experiment, we calculated the error rate in the whole experiment and then applied it to every male to avoid discarding genuine haplotypes that belong to only one or few males. To run the script, generate other "lista.txt" file with the path of all the "*_hap.fas" files by executing "$ ls *hap.fas > lista.txt".
```
$ python 04umbralizator.py
```
The script will ask the number of threads to execute MAFFT. For each sample, it generate two aligned FASTA files fot the ITS2 and the contiguous genes named "*inner_hap.fas" and "*outer_hap.fas", respectively. In addition, you get a file with the frequency of the more abundant haplotype ("umbrales.txt") and a FASTA file with all the haplotypes for all the found ITS2 and genes haplotypes ("hap_total.fas").
5) The reads were then classified according to the ITS2 haplotype to which they belong.
```
$ python 05.haplotyper.py
```

# References

* Bragg, L., Stone, G., Imelfort, M., Hugenholtz, P. and Tyson, G.W. (2012). Fast, accurate error-correction of amplicon pyrosequences using Acacia. Nature Methods 9: 425-426.
* Edgar, R.C., Haas, B.J., Clemente, J.C., Quince, C. and Knight, R. (2011). UCHIME improves sensitivity and speed of chimera detection. Bioinformatics 27(16): 2194-2200.
* Faircloth, B.C. and Glenn, T.C. (2012). Not all sequence tags are created equal: designing and validating sequence identification tags robust to indels. PloS One 7:e42543. doi: 10.1371/journal.pone.0042543.
* Katoh, K. and Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol 30: 772-780.
* Koetschan, C., Förster, F., Keller, A., Schleicher, T., Ruderisch, B. and Schwarz, R. (2010) The ITS2 Database III - Sequences and structures for phylogeny. Nucleic Acids Res 38(Database issue): D275-D279.
* Rice, P., Longden, I. and Bleasby, A. (2000). EMBOSS: the European molecular biology open software suite. Trends in Genetics 16: 276-277.
* Teruel, M., Ruiz-Ruano, F.J., Marchal, J.A., Sánchez, A., Cabrero, J., Camacho, J.P.M. et al. (2014). Disparate molecular evolution of two types of repetitive DNAs in the genome of the grasshopper Eyprepocnemis plorans. Heredity 112:531-42.
