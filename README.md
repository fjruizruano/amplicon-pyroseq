amplicon-pyroseq
================

Processing and analyzing pyrosequencing reads from amplicons

# Dependencies

* EMBOSS
* Acacia
* Uchime
* MAFFT

# Pipeline


We wrote a series of custom Python scripts and used them for counting the number of reads corresponding to each sequence (haplotype) in each sample. This was done in five consecutive steps:

1) To select the ITS2 amplicons, local alignments of the initial half of the sequencing reads against the forward and reverse primers were carried out using the Smith-Waterman algorithm, as implemented in the EMBOSS suite (Rice et al., 2000). We only considered reads with high identity with one of the primers. Additionally, this alignment served to determine the orientation of each read in respect of the aligned primer, since we had reads in both orientations because the 454 adaptors were not included in the primers. We substituted the reads matched to the reverse primer with their reverse-complementary in order to have all the amplicons in the same orientation.
```
$ python 01split_genes.py
```
2) We assigned each read to a given sample according to its tag combination using local alignments of the complete reads againts the tagged primers. 6mer Tags were designed using using EDITTAG (Faircloth and Glenn, 2012), differing in four or more nucleotides. To avoid missassignation, we sorted the alignments from higher to lower similarity and assigned the sample to a given read only if the best alignment showed less than four sequence differences, (i.e. the minimal number of differences between tags) in both tagged primers. Selected reads were stored in different files according to their combination of tags. The alignments also allowed us to trim tags and primers from the sequences.
```
$ python 02split_tag.py
```
3) We used the Acacia software (Bragg et al., 2012) to correct typical 454 pyrosequencing errors, i.e. substitutions and indels in homopolymeric regions. We then generated a file with all found haplotypes per sample and their frequency, and we searched for chimeric sequences with UCHIME (Edgar et al., 2011) using the default options for the de novo algorithm.
```
$ python 03clean_amp.py
```
4) To annotate the 5.8S, ITS2 and 28S regions in all the haplotype files, per sample, we aligned the sequences using MAFFT v7 (Katoh and Standley, 2013) with LINSI options including an additional E. plorans ITS2 sequence (accession number JN811835.1, corresponding with the ref.fas file) as a reference, annotating its ITS2 region with the ITS2 Database III annotation tool (Koetschan et al., 2010).
```
$ python 04umbralizator.py
```
5) Since the 454 reads included partial regions of the 5.8S and 28S rRNA genes, summing up 123 nt, in addition to the ITS2 sequence, we used these partial coding regions as internal control for sequencing errors. This provided a way to avoid false positives in identifying genuine ITS2 haplotypes (i.e. variants for the ITS2 region sequence) in the gDNA. For stringency, we considered as sequencing errors all the variation found in these 123 nt, compared to the sequence found by Teruel et al. (2014) which is conserved in all E. plorans rDNA sequences (accession number JN811835.1). The proportion of reads carrying any variation in respect to the conserved coding sequence was thus our estimate of the maximum error rate of the experiment. In order to select the genuine ITS2 haplotypes in the different samples of the experiment, we calculated the error rate in the whole experiment and then applied it to every male to avoid discarding genuine haplotypes that belong to only one or few males. The reads were then classified according to the ITS2 haplotype to which they belong.
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
