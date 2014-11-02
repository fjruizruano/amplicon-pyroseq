#!/usr/bin/python
from Bio import SeqIO, AlignIO
from subprocess import call

threads = raw_input("number of threads: ")

files = open("lista.txt").readlines()
ref = open("ref.fas").readlines()
ref = ref[1][:-1]
ref = ref.split("-")

print ref

for file in files:
	name = file[:-5]
	call("cat ref.fas %s > %s" % (file[:-1], name+"_ref.fas"), shell=True)
	call("mafft-linsi --thread %s %s > %s" % (threads, name+"_ref.fas", name+"_ref_m.fas"), shell=True)
	ali = AlignIO.read(open(name+"_ref_m.fas"), "fasta")
	print ali[0]
	ali_ref = str(ali[0].seq)
	ali_ref = ali_ref.upper()

	i = -1
	
	for col in range(0,len(ali_ref)):
		print len(ali_ref)
		if ali_ref[col] != "-":
			i += 1
			if i == 0:
				begin1 = col

			elif i == len(ref[0])-1:
				end1 = col

			elif i == len(ref[0]):
				begin2 = col

			elif i == len(ref[0])+len(ref[1])-1:
				end2 = col

	subali_outer = ali[1:,begin1:end1+1] + ali[1:,begin2:end2+1]
	subali_inner = ali[1:,end1+1:begin2]

	SeqIO.write(subali_outer, name+"_outer.fas", "fasta")
	SeqIO.write(subali_inner, name+"_inner.fas", "fasta")
