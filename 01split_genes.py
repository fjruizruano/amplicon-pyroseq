from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import operator

class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
        
    def __iter__(self):
        return self
    
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
        
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

#its1 = "TAGGTGAACCTGCGGAAG"
#its2 = "TCCACCGTCCTGGGTGATCTT"
its3 = "GTCGATGAAGAACGCAGC"
its4 = "ATATGCTTAAATTCAGCGGG"
hist3 = "TCGTACATGTGTTATGTGAGC" #scarf
hist4 = "TGCTGCCGGCGACAGGTTCACGTA"  #scarr
#satL = "CGCATTTCTGCCGCCTG"
#satR = Seq("GCACTGCTTTCCAGATAT")

#primers = [its1, its2, its3, its4, satL, satR]
#primers_name = ["its1", "its2", "its3", "its4", "satL", "satR"]

primers = [its3, its4, hist3, hist4 ]
primers_name = ["its3", "its4", "hist3", "hist4"]

secu = ParseFastQ("pcr5.fq")

#out12 = open("its12.fq", "w")
out34 = open("its34.fq", "w")
outhi = open("hist34.fq", "w")
#outLR = open("satLR.fq", "w")
outNM = open("noMat.fq", "w")

j = 0

for a in secu:
	j += 1

total = j
print total
total = total * 1.0

secu = ParseFastQ("pcr5.fq")

#uno = 0
dos = 0
his = 0
#sat = 0
nma = 0

for line in secu:

	j -= 1

	minu = list(line[1])[6:40] # tag IMPORTANTE REVISAR

	minu = "".join(minu)
	mayu = minu
#	mayu = str.upper(minu)
	print "MAYUS: " + mayu

	scores = {}
	for primer_name in primers_name:
		scores[primer_name] = 0

	count_name = 0

	for primer in primers:
		x = pairwise2.align.localms(primer, mayu, 1, -1, -1, -.1) 
		x = x[0][2]
		x = int(x)
		scores[primers_name[count_name]] = x
		count_name += 1
	sort_scores = sorted(scores.iteritems(), key=operator.itemgetter(1))
	print sort_scores
	best_score = sort_scores[3] # IMPORTANTE REVISAR
	print str(round((j/total*100), 2))

	if best_score[1] > 14: #IMPORTANTE REVISAR

		if best_score[0] == "its3":
			dos += 1
			out34.write("%s_its34_F_seq%s\n%s\n%s\n%s\n" % (line[0], str(dos),line[1], line[2], line[3]))

		elif best_score[0] == "its4":
			dos += 1
			line_rc = Seq(line[1], generic_dna)
			line_rc = str(line_rc.reverse_complement())
			line_q = list(line[3])
			line_q.reverse()
			line_q = "".join(line_q)
			out34.write("%s_its34_R_seq%s\n%s\n%s\n%s\n" % (line[0], str(dos), line_rc, line[2], line_q))

		elif best_score[0] == "hist3":
			his += 1
			outhi.write("%s_hist34_F_seq%s\n%s\n%s\n%s\n" % (line[0], str(his),line[1], line[2], line[3]))

		elif best_score[0] == "hist4":
			his += 1
			line_rc = Seq(line[1], generic_dna)
			line_rc = str(line_rc.reverse_complement())
			line_q = list(line[3])
			line_q.reverse()
			line_q = "".join(line_q)
			outhi.write("%s_hist34_R_seq%s\n%s\n%s\n%s\n" % (line[0], str(his), line_rc, line[2], line_q))

	else:
		nma += 1
		outNM.write("%s_noMat_seq%s\n%s\n%s\n%s\n" % (line[0], str(nma),line[1], line[2], line[3]))

#out12.close()
out34.close()
outhi.close()
#outLR.close()
outNM.close()
