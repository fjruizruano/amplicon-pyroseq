from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
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

its1_tag2 = "AGGCTATAGGTGAACCTGCGGAAG"
its1_tag3 = "CCATTGTAGGTGAACCTGCGGAAG"
its1_tag4 = "CTCGAATAGGTGAACCTGCGGAAG"
its1_tag5 = "TGTACGTAGGTGAACCTGCGGAAG"

its2_tag2 = Seq("AGGCTATCCACCGTCCTGGGTGATCTT")
its2_tag2 = str(its2_tag2.reverse_complement())
its2_tag3 = Seq("CCATTGTCCACCGTCCTGGGTGATCTT")
its2_tag3 = str(its2_tag3.reverse_complement())
its2_tag4 = Seq("CTCGAATCCACCGTCCTGGGTGATCTT")
its2_tag4 = str(its2_tag4.reverse_complement())
its2_tag5 = Seq("TGTACGTCCACCGTCCTGGGTGATCTT")
its2_tag5 = str(its2_tag5.reverse_complement())


its3_tag1 = "AACAGCGTCGATGAAGAACGCAGC"
its3_tag2 = "AGGCTAGTCGATGAAGAACGCAGC"
its3_tag3 = "CCATTGGTCGATGAAGAACGCAGC"
its3_tag4 = "CTCGAAGTCGATGAAGAACGCAGC"
its3_tag5 = "TGTACGGTCGATGAAGAACGCAGC"
its3_tag6 = "AGAATAGTCGATGAAGAACGCAGC"

its4_tag1 = Seq("AACAGCATATGCTTAAATTCAGCGGG")
its4_tag1 = str(its4_tag1.reverse_complement())
its4_tag2 = Seq("AGGCTAATATGCTTAAATTCAGCGGG")
its4_tag2 = str(its4_tag2.reverse_complement())
its4_tag3 = Seq("CCATTGATATGCTTAAATTCAGCGGG")
its4_tag3 = str(its4_tag3.reverse_complement())
its4_tag4 = Seq("CTCGAAATATGCTTAAATTCAGCGGG")
its4_tag4 = str(its4_tag4.reverse_complement())
its4_tag5 = Seq("TGTACGATATGCTTAAATTCAGCGGG")
its4_tag5 = str(its4_tag5.reverse_complement())
its4_tag6 = Seq("AGAATAATATGCTTAAATTCAGCGGG")
its4_tag6 = str(its4_tag6.reverse_complement())


hist3_tag1 = "AACAGCCTTGGTGGCGAGCTGTTTGCG"
hist3_tag2 = "AGGCTACTTGGTGGCGAGCTGTTTGCG"
hist3_tag3 = "CCATTGCTTGGTGGCGAGCTGTTTGCG"
hist3_tag4 = "CTCGAACTTGGTGGCGAGCTGTTTGCG"
hist3_tag5 = "TGTACGCTTGGTGGCGAGCTGTTTGCG"

hist4_tag1 = Seq("AACAGCCCCTTGCCGAGCCCCTTTCC")
hist4_tag1 = str(hist4_tag1.reverse_complement())
hist4_tag2 = Seq("AGGCTACCCTTGCCGAGCCCCTTTCC")
hist4_tag2 = str(hist4_tag2.reverse_complement())
hist4_tag3 = Seq("CCATTGCCCTTGCCGAGCCCCTTTCC")
hist4_tag3 = str(hist4_tag3.reverse_complement())
hist4_tag4 = Seq("CTCGAACCCTTGCCGAGCCCCTTTCC")
hist4_tag4 = str(hist4_tag4.reverse_complement())
hist4_tag5 = Seq("TGTACGCCCTTGCCGAGCCCCTTTCC")
hist4_tag5 = str(hist4_tag5.reverse_complement())


scarf_tag4 = "CTCGAATCGTACATGTGTTATGTGAGC"
scarf_tag3 = "CCATTGTCGTACATGTGTTATGTGAGC"
scarf_tag2 = "AGGCTATCGTACATGTGTTATGTGAGC"

scarr_tag4 = Seq("CTCGAATGCTGCCGGCGACAGGTTCACGTA")
scarr_tag4 = str(scarr_tag4.reverse_complement())
scarr_tag3 = Seq("CCATTGTGCTGCCGGCGACAGGTTCACGTA")
scarr_tag3 = str(scarr_tag3.reverse_complement())
scarr_tag2 = Seq("AGGCTATGCTGCCGGCGACAGGTTCACGTA")
scarr_tag2 = str(scarr_tag2.reverse_complement())


#satL_tag1 = "AACAGCCGCATTTCTGCCGCCTGTGGCGCTACATT"
#satL_tag3 = "CCATTGCGCATTTCTGCCGCCTGTGGCGCTACATT"
#satL_tag4 = "CTCGAACGCATTTCTGCCGCCTGTGGCGCTACATT"
#satL_tag5 = "TGTACGCGCATTTCTGCCGCCTGTGGCGCTACATT"

#satR_tag1 = Seq("AACAGCGCACTGCTTTCCAGATATCACACTAAAATG")
#satR_tag1 = str(satR_tag1.reverse_complement())
#satR_tag2 = Seq("AGGCTAGCACTGCTTTCCAGATATCACACTAAAATG")
#satR_tag2 = str(satR_tag2.reverse_complement())
#satR_tag3 = Seq("CCATTGGCACTGCTTTCCAGATATCACACTAAAATG")
#satR_tag3 = str(satR_tag3.reverse_complement())
#satR_tag5 = Seq("TGTACGGCACTGCTTTCCAGATATCACACTAAAATG")
#satR_tag5 = str(satR_tag5.reverse_complement())

key = raw_input("Introduce gene name: ")
dist = raw_input("Introduce distance: ")
dist = int(dist)

# lista de primers
if key == "its1":
	forwards = [its1_tag2, its1_tag3, its1_tag4, its1_tag5]
	forwards_names = ["its1_tag2", "its1_tag3", "its1_tag4", "its1_tag5"]
	len_for = len(forwards_names)
	reverses = [its2_tag2, its2_tag3, its2_tag4, its2_tag5]
	reverses_names = ["its2_tag2", "its2_tag3", "its2_tag4", "its2_tag5"]
	len_rev = len(reverses_names)
	input = "its12.fq"

elif key == "its2":
	forwards = [its3_tag1, its3_tag2, its3_tag3, its3_tag4, its3_tag5, its3_tag6]
	forwards_names = ["its3_tag1", "its3_tag2", "its3_tag3" , "its3_tag4", "its3_tag5", "its3_tag6"]
	len_for = len(forwards_names)
	reverses = [its4_tag1, its4_tag2, its4_tag3, its4_tag4, its4_tag5, its4_tag6]
	reverses_names = ["its4_tag1", "its4_tag2", "its4_tag3", "its4_tag4", "its4_tag5", "its4_tag6"]
	len_rev = len(reverses_names)
	input = "its34.fq"

elif key == "hist":
	forwards = [hist3_tag1, hist3_tag2, hist3_tag3, hist3_tag4, hist3_tag5]
	forwards_names = ["hist3_tag1", "hist3_tag2", "hist3_tag3", "hist3_tag4", "hist3_tag5"]
	len_for = len(forwards_names)
	reverses = [hist4_tag1, hist4_tag2, hist4_tag3, hist4_tag4, hist4_tag5]
	reverses_names = ["hist4_tag1", "hist4_tag2", "hist4_tag3", "hist4_tag4", "hist4_tag5"]
	len_rev = len(reverses_names)
	input = "hist34.fq"

elif key == "scar":
	forwards = [scarf_tag4, scarf_tag3, scarf_tag2]
	forwards_names = ["scarf_tag4", "scarf_tag3", "scarf_tag2"]
	len_for = len(forwards_names)
	reverses = [scarr_tag4, scarr_tag3, scarr_tag2]
	reverses_names = ["scarr_tag4", "scarr_tag3", "scarr_tag2"]
	len_rev = len(reverses_names)
	input = "scarfr.fq"

#elif key == "sat":
#	forwards = [satL_tag1, satL_tag3, satL_tag4, satL_tag5]
#	forwards_names = ["satL_tag1", "satL_tag3", "satL_tag4", "satL_tag5"]
#	reverses = [satR_tag1, satR_tag2, satR_tag3, satR_tag5]
#	reverses_names = ["satR_tag1", "satR_tag2", "satR_tag3", "satR_tag5"]
#	input = "satLR.fq"

secu = ParseFastQ(input)

j = 0

for a in secu:
	j += 1

total = j
print total
total = total * 1.0

secu = ParseFastQ(input)

i = 0

currout = None
data = None
seqname2file = dict()
counters = dict()

for line in secu:

	j -= 1

	minu = line[1] 
	mayu = str.upper(minu)
	qual = line[3]

	# primer forward
	scores_f = {}
	for forward_name in forwards_names:
		scores_f[forward_name] = 0
	count_name_f = 0
	for forward in forwards:
		# Alineamiento local del primer forward
		x = pairwise2.align.localms(forward, mayu[0:len(mayu)/2], 1, -1, -1, -.1) 
		xs = x[0][2] # puntuacion
		xb = x[0][3] # inicio
		xl = x[0][4] # fin
		xlist = [xs, xb, xl]
		scores_f[forwards_names[count_name_f]] = xlist
		count_name_f +=1
	sort_scores_f = sorted(scores_f.iteritems(), key=operator.itemgetter(1))
	bsf = sort_scores_f[len_for-1] # best score for forward primer

	# primer reverse
	scores_r = {}
	for reverse_name in reverses_names:
		scores_r[reverse_name] = 0
	count_name_r = 0
	for reverse in reverses:
		y = pairwise2.align.localms(reverse, mayu[len(mayu)/2:-1], 1, -1, -1, -.1)
		ys = y[0][2]
		yb = y[0][3]
		yl = y[0][4]
		ylist = [ys, yb, yl]
		scores_r[reverses_names[count_name_r]] = ylist
		count_name_r += 1
	sort_scores_r = sorted(scores_r.iteritems(), key=operator.itemgetter(1))
	bsr = sort_scores_r[len_rev-1] # best score for reverse primer


#	print str(bsf[1][0]) + ", " + str(bsr[1][0]) + " | " + str(round((j/total*100), 2))

	# Si hay 4 cambios entre los dos primers, lo coge eliminando los primers
	if bsf[1][0] >= len(forward)-dist and bsr[1][0] >= len(reverse)-dist:
		
		tag = str(bsf[0] + "-" + bsr[0])

		if tag not in seqname2file:
			filename = "out_%s.fq" % tag
			seqname2file[tag] = open(filename, "w")

		if tag not in counters:
			counters[tag] = 1
		else:
			counters[tag] += 1
		data = "%s_seq_%s\n%s\n%s\n%s\n" % (line[0], str(counters[tag]), mayu[bsf[1][2]:(len(mayu)/2)+bsr[1][1]], line[2], qual[bsf[1][2]:(len(mayu)/2)+bsr[1][1]])
		currout = seqname2file[tag]
		currout.write(data)
#		print counters[tag]
#		print data	
#		print counters

for f in seqname2file.values():
	f.close()
