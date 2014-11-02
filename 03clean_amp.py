from Bio import SeqIO
import operator
import os
from subprocess import call

if not os.path.exists("acacia_files"):
	os.makedirs("acacia_files")

files = open("lista.txt", "r").readlines()
for file in files:
	acacia_conf = open("acacia.conf", "r").readlines()
	acacia_conf[6] = "FASTQ_LOCATION=./" + file
	acacia_conf[18]= "OUTPUT_PREFIX="  + file[:-4] + "\n"
	acacia_conf_w = open("acacia2.conf", "w")
	acacia_conf_w.write("".join(acacia_conf))
	acacia_conf_w.close()
	call("acacia_1.52b0 -c ./acacia2.conf", shell=True)
	secu = SeqIO.parse(open("./acacia_files/%s_all_tags.seqOut" % file[:-4]), "fasta")

	dictio = {}
	i = 0

	for seque in secu:
		look = str(seque.seq) in dictio
		if look == True:
			dictio[str(seque.seq)] = dictio[str(seque.seq)] + 1
		elif look == False:
			dictio[str(seque.seq)] = 1
			i += 1

	dictio = sorted(dictio.iteritems(), key=operator.itemgetter(1))
	dictio.reverse()

	output_handle = open(file[:-4] + "_hap.fas", "w")

	i = 0

	for sets in dictio:
		i += 1
		output_handle.write(">Seq_%s;size=%s\n%s\n" % (str(i), str(sets[1]), sets[0]))
		
	print dictio
	print i

	call("usearch6.0.307_i86linux32 -uchime_denovo %s -uchimeout %s" % (file[:-4]+"_hap.fas", file[:-4]+"_hap_chime.out" ), shell=True)

	# CONTINUAR ELIMINANDO CHIMERAS
