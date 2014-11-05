#!/usr/bin/python
#*-* coding: utf-8 *-*

from Bio import SeqIO

files_chr = open("lista.txt").readlines() #cargar lista ficheros fasta 
file_hap = SeqIO.parse(open("hap_output.fas"), "fasta") #cargar fichero haplotipos

count_dict = {}
count_table = {}
header = ""

#i = 0

#crear diccionarios con los haplotipos
for sequen in file_hap:
#	i +=1
	count_dict[str(sequen.seq)] = 0
	count_table[str(sequen.seq)] = [sequen.id]

i = 0

#para cada fichero de secuencias
for line in files_chr:
	cut = line.rfind("/")
	header += "\t%s" % line[cut+1:-5]

	i += 1

	chr = SeqIO.parse(open("%s" % line[:-1]), "fasta") #fichero

	#para cada secuencia en el fichero 
	for secu in chr:
		#ver si la secuencia está en el diccionario, sumamos su size
		look_seq = str(secu.seq) in count_dict
		if look_seq == True:
			x = secu.id
			x = x.split("size=")
			try:
				counter = int(x[1])
				count_dict[str(secu.seq)] += counter
			except:
				print "error en secuencia %s fichero %s" % (str(secu.id), str(i))

	#añadimos al diccionario los contemos y lo reinicamos los conteos parciales
	for a in count_dict:
		count_table[a].append(count_dict[a])
		count_dict[a] = 0

#grabamos todos los resultados
out = open("output.txt", "w")
out.write(header + "\n")

for a in count_table:
	print count_table[a]
	out.write("\t".join(map(str,count_table[a]))+"\n")
out.close()

# crea lista de muestras
desample = header.split("\t")
desample = desample[1:]

# crea diccionario de muestras
dict_seq = {}
for sample in desample:
	dict_seq[sample] = []

hap_ali = SeqIO.parse(open("haplotipos_individuos_m.fas"), "fasta")

dict_hap = {}
for hap in hap_ali:
	dict_hap[hap.id] = str(hap.seq)

# rellena diccionario de muestras con todas las secuencias
for a in count_table:
	counts = count_table[a]
	haplo = counts[0] # id del haplotipo
	counts = counts[1:] # conteo del hap para cada muestra
	a = dict_hap[haplo] # haplotipo con gap
	i = 0
	for count in counts:
		add = [a]*count
		sample = desample[i]
		dict_seq[sample] += add
		i += 1
#print counts

k = -1
for x in dict_seq:
	k += 1
	l = 0
	out_hap = open(x+".fas", "w")
	for seq in dict_seq[x]: 
		l += 1
		out_hap.write(">%s_%s\n%s\n" % (str(x), str(l), str(seq)))
	out_hap.close()
