#!/usr/bin/python
#*-* coding: utf-8 *-*

from Bio import SeqIO, AlignIO
from subprocess import call
import operator

threads = raw_input("number of threads: ")

files = open("lista.txt").readlines()
ref = open("ref.fas").readlines()
ref = ref[1][:-1]
ref = ref.split("-")

haplotypes = []
#dict_umbral = {}
dict_total = {}

out_umbral = open("umbrales.txt", "w")

for file in files:
	name = file[:-5]
	call("cat ref.fas %s > %s" % (file[:-1], name+"_ref.fas"), shell=True)
	call("mafft-linsi --thread %s %s > %s" % (threads, name+"_ref.fas", name+"_ref_m.fas"), shell=True)
	ali = AlignIO.read(open(name+"_ref_m.fas"), "fasta")
	ali_ref = str(ali[0].seq)
	ali_ref = ali_ref.upper()

	i = -1
	
	for col in range(0,len(ali_ref)):
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

	trimmed = [name+"_outer.fas", name+"_inner.fas"]

	SeqIO.write(subali_outer, trimmed[0], "fasta")
	SeqIO.write(subali_inner, trimmed[1], "fasta")

	#reducir a haplotipos

	number = 0
	for trim in trimmed:

		number += 1
		file = trim
#		print file
	
		#abrimos el fichero
		secu = SeqIO.parse(file, "fasta")

		#contamos el numero total de secuencias por fichero
		i = 0
		dict_ind = {}
		for sequen in secu:
			size = sequen.id
			size = size.split("size=")
			size = int(size[-1])
#			if size > 1: #solo si la sec esta mas de una vez
#				i += size
	
			#añadimos secuencia y abundancia
			sec = str(sequen.seq)
			sec = sec.replace("-", "")
			look = sec in dict_ind
			if look == True:
				dict_ind[sec] += size
			elif look == False:
				dict_ind[sec] = size
			
			look2= sec in dict_total
			if look2== True:
				dict_total[sec] += size
			elif look2== False:
				dict_total[sec] = size
	
		#diccionario de umbrales por fichero
#		dict_umbral[file] = i

		#ordena de mayor a menor abundancia por fichero
		dict_ind = sorted(dict_ind.iteritems(), key=operator.itemgetter(1))
		dict_ind.reverse()
	
		#para cada fichero añadimos un haplotipo hasta el umbral de 95%
#		sum = 0
		for hap in dict_ind:
#			if sum < 0.9999*i:
				look = hap[0] in haplotypes
				if look == False:
					haplotypes.append(hap[0])
#				sum += hap[1]
	
		#escribe ficheros de haplotipos nuevos
		file_ind = open("%s_hap.fas" % file[:-4], "w")
		h = 0
		total = 0
		first = 0
		for haplo in dict_ind:
			h += 1
			total += haplo[1]
			if h == 1:
				first = haplo[1]
			file_ind.write(">Hap_%s;size=%s\n%s\n" % (str(h), str(haplo[1]), haplo[0]))
		file_ind.close()

		if number%2 == 1:
			out_umbral.write(file + "\t" + str(1.0*first/total) + "\n")
	
#escribe fichero de todos los haplotipos seleccionados
#out = open("hap_output.fas", "w")
#j = 0
#for seq in haplotypes:
#	j += 1
#	out.write(">Hap_%s\n%s\n" % (str(j), seq))
#out.close()


dict_total = sorted(dict_total.iteritems(), key=operator.itemgetter(1))
dict_total.reverse()

out = open("hap_total.fas", "w")
h=0
for haplo in dict_total:
	h += 1
	out.write(">Hap%s_size=%s\n%s\n" % (str(h), str(haplo[1]), haplo[0]))
out.close()
out_umbral.close()
