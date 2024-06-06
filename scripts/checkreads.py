from skbio import alignment,sequence
from Bio import SeqIO,SeqRecord,Seq
from skbio import alignment
from statistics import mean,stdev

import sys
def reverse_complement(dna):
    """

    :param dna:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "-": "-", "N": "N"}
    return ''.join([complement[base] for base in dna[::-1]])

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_sw_match(primer,read):
    readseq = sequence.DNA(read)
    primerseq = sequence.DNA(primer)
    onegapor2missmatch = (len(primer) * 2) - 11
    res = alignment.local_pairwise_align_ssw(readseq, primerseq, gap_open_penalty=5, gap_extend_penalty=2, match_score=2,
                                             mismatch_score=-3,score_filter=onegapor2missmatch-1)


    if res:
        return res[2]
    return False

def read_cleanup(r,adapters):
    for ad in adapters:
        adseq = str(adapters[ad].seq)
        rseq = str(r.seq)
        if adseq in rseq:
            adpos = rseq.find(adseq)
            newrseq = rseq[:adpos]
            return newrseq

def check_for_ag(r):
    if "AAAAAAGGGGGGGGGG" in str(r.seq):
        return True
    else:
        return False

def remove_polyg(r):
    polyg = "GGGGG"
    r1seq = str(r.seq)
    if polyg in r1seq:
        polygpos = r1seq.find(polyg)
        newr1seq = r1seq[:polygpos]
        r = SeqRecord.SeqRecord(Seq.Seq(newr1seq), id=r.id, description=r.description,
                                      name=r.name, letter_annotations={
                "phred_quality": r.letter_annotations["phred_quality"][:len(newr1seq)]})
    return r

def remove_adseq(r,adseq):
    rseq = str(r.seq)
    revadseq = reverse_complement(adseq)
    found = False
    if adseq in rseq:
        adpos = rseq.find(adseq)
        newrseq = rseq[:adpos]
        r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                      name=r.name, letter_annotations={
                "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})
        found = True
    if not found:
        if revadseq in rseq:
            adpos = rseq.find(revadseq)
            newrseq = rseq[:adpos]
            r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                          name=r.name, letter_annotations={
                    "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})
            found = True
    if not found:
        len_ad = len(adseq)
        hdist = hamming_distance(adseq, rseq[-1*len_ad:])
        if hdist < 10:
            newrseq = rseq[:len_ad]
            r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                    name=r.name, letter_annotations={
                    "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})
            found = True

    if not found:
        len_ad = len(adseq)
        hdist = hamming_distance(revadseq, rseq[-1*len_ad:])
        if hdist < 10:
            newrseq = rseq[:len_ad]
            r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                    name=r.name, letter_annotations={
                    "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})
            found = True

    return r

class readpair():
    def __init__(self):
        self.r1 = False
        self.r2 = False
        self.olread = False
    def add1(self,name,r1):
        self.name = name
        self.r1 = r1
        self.r1len = len(r1.seq)
        self.r1hasag = check_for_ag(self.r1)
    def add2(self,name,r2):
        self.name = name
        self.r2 = r2
        self.r2len = len(r2.seq)
        self.r2hasag = check_for_ag(self.r2)
    def add_ol(self,olread):
        self.olread = olread
    def readcleanup(self,adapters):

        for ad in adapters:
            adseq = str(adapters[ad].seq)
            self.r1 = remove_adseq(self.r1,adseq)
            self.r1len = len(self.r1.seq)
            self.r2 = remove_adseq(self.r2,adseq)
            self.r2len = len(self.r2.seq)

    def add_primer_type(self,p1,p2):
        self.p1 = p1
        self.p2 = p2
        if self.p1 and self.p2:
            self.mtype = (self.p1.id, self.p2.id)
            if self.p1.id.startswith("rc_"):
                p1gene = self.p1.id.replace("rc_","")
            else:
                p1gene = self.p1.id
            p1gene = p1gene.split("_")[0]
            if self.p2.id.startswith("rc_"):
                p2gene = self.p2.id.replace("rc_", "")
            else:
                p2gene = self.p2.id
            p2gene = p2gene.split("_")[0]
            if p1gene == p2gene and self.p1.id.endswith("_f") and self.p2.id.endswith("_r") and "rc" not in self.p1.id and "rc" not in self.p2.id:
                self.type = "correct"
            elif p1gene == p2gene and self.p1.id.endswith("_f") and self.p2.id.endswith("_f"):
                self.type = "possible_homodimer"
            elif p1gene == p2gene and self.p1.id.endswith("_r") and self.p2.id.endswith("_r"):
                self.type = "possible_homodimer"
            elif p1gene == p2gene and self.p1.id.endswith("_r") and self.p2.id.endswith("_f"):
                self.type = "inverted"
            else:
                self.type = "mixed"

        elif self.p1 or self.p2:
            self.type = "single"
            if self.p1:
                self.mtype = (self.p1.id, False)
            else:
                self.mtype = (False, self.p2.id)
        else:
            self.type = "none"
            self.mtype = (False, False)
    def check_products(self,nonprimertargets,targets,snpdata):
        self.r1prodmatch = False
        self.r2prodmatch = False
        self.olmatch = False
        self.snptype = False
        if self.p1 and self.p2 and self.olread:
            if self.p1.id.startswith("rc_"):
                p1gene = self.p1.id.replace("rc_","")
            else:
                p1gene = self.p1.id
            gene = p1gene.split("_")[0]
            gene_expected_prod = targets[gene + "_amplicon"]
            hdist = hamming_distance(str(self.olread.seq), str(gene_expected_prod.seq))
            lenratio = len(str(self.olread.seq))/len(str(gene_expected_prod.seq))

            if float(hdist) < float(len(gene_expected_prod.seq)/10) and lenratio > 0.9:
                self.olmatch = True
                snppos = int(snpdata[gene][0])
                snpref = snpdata[gene][1]
                snpalt = snpdata[gene][2]
                olseq = str(self.olread.seq)

                try:
                    if olseq[snppos-1] == snpref:
                        self.snptype = (gene,"ref",snppos,snpref)
                    elif olseq[snppos-1] == snpalt:
                        self.snptype = (gene,"alt",snppos,snpalt)
                    else:
                        self.snptype = (gene,"other",snppos,olseq[snppos-1])
                except Exception as e:
                    print("\n\nERROR\n")
                    print(e)
                    print("\n\n")
        if self.p1:
            if self.p1.id.startswith("rc_"):
                p1gene = self.p1.id.replace("rc_","")
            else:
                p1gene = self.p1.id
            p1gene = p1gene.split("_")[0]
            p1len = len(self.p1.seq)
            p1gene_expected_prod = nonprimertargets[p1gene+"_noprimers"]
            if len(self.r1.seq) > (len(p1gene_expected_prod)+p1len):
                noprimer_read_seq = str(self.r1.seq)[p1len:len(p1gene_expected_prod)]
            else:
                noprimer_read_seq = str(self.r1.seq)[p1len:]
            noprimer_p1gene_expected_prod = str(p1gene_expected_prod.seq)[:len(noprimer_read_seq)]
            hdist = hamming_distance(noprimer_read_seq,noprimer_p1gene_expected_prod)
            if hdist < 10 and len(noprimer_read_seq) > 40:
                self.r1prodmatch = True

        if self.p2:
            if self.p2.id.startswith("rc_"):
                p2gene = self.p2.id.replace("rc_", "")
            else:
                p2gene = self.p2.id
            p2gene = p2gene.split("_")[0]
            p2len = len(self.p2.seq)
            p2gene_expected_prod = nonprimertargets[p2gene+"_noprimers"].reverse_complement()
            if len(self.r2.seq) > (len(p2gene_expected_prod)+p2len):
                noprimer_read_seq = str(self.r2.seq)[p2len:len(p2gene_expected_prod)]
            else:
                noprimer_read_seq = str(self.r2.seq)[p2len:]

            noprimer_p2gene_expected_prod = str(p2gene_expected_prod.seq)[:len(noprimer_read_seq)]

            hdist = hamming_distance(noprimer_read_seq,noprimer_p2gene_expected_prod)

            if hdist < 10 and len(noprimer_read_seq) > 40:
                self.r2prodmatch = True
        if self.type == "correct" and (len(self.r1.seq) < 50 or len(self.r2.seq) < 50):
            self.type = "short_reads"
        elif self.olread:
            if not self.olmatch and self.type == "correct":
                self.type = "ovelap_wrong"
                
    def check_primer_dimer(self):
        self.primerdimer = []
        if self.r1len == self.r2len:
            if self.p1 and self.p2:
                p1len = len(self.p1.seq)
                p2len = len(self.p2.seq)
                if self.r1len < (p1len+p2len):
                    self.primerdimer.append("dimer_readlengths")
                elif self.r1len < 130:
                    self.primerdimer.append("Short_readlength")
        if self.olread:
            ollen = len(self.olread.seq)
            if self.p1 and self.p2:
                p1len = len(self.p1.seq)
                p2len = len(self.p2.seq)
                if ollen < (p1len + p2len):
                    self.primerdimer.append("dimer_overlap_length")
                elif ollen < 130:
                    self.primerdimer.append("Short_overlaplength")
            elif self.p1:
                p1len = len(self.p1.seq)
                if ollen <= p1len:
                    self.primerdimer.append(self.p1.id + "_only")
            elif self.p2:
                p2len = len(self.p2.seq)
                if ollen <= p2len:
                    self.primerdimer.append(self.p2.id + "_only")

    def to_out_tab(self):
        concatm = str(self.mtype[0])+":"+str(self.mtype[1])
        if self.p1:
            p1=str(self.p1.id)
        else:
            p1 = str(self.p1)
        if self.p2:
            p2 = str(self.p2.id)
        else:
            p2 = str(self.p2)
        if self.olread:
            ollen = len(self.olread.seq)
        else:
            ollen = ""
        return "{}\n".format("\t".join([self.name,self.type,p1,p2,concatm,str(self.r1hasag),str(self.r2hasag),str(self.r1prodmatch),str(self.r2prodmatch),str(self.olmatch),str(self.r1len),str(self.r2len),str(ollen),str(":".join(self.primerdimer))]))
    def __str__(self):
        return f"{self.name}\t{self.type}\t{self.mtype}\t{self.orient}"


def checkreads(r1,r2,olreads,primers,nonprimertargets,targets,adapters,output,snpdata):

    read_details = output.replace(".txt","_readinfo.txt")
    snptyping = output.replace(".txt", "_snps.txt")
    inread = SeqIO.parse(r1,"fastq")
    inread2 = SeqIO.parse(r2,"fastq")
    inoverlaps = SeqIO.parse(olreads,"fastq")


    primers = SeqIO.to_dict(SeqIO.parse(primers,"fasta"))
    nonprimertargets = SeqIO.to_dict(SeqIO.parse(nonprimertargets,"fasta"))
    targets = SeqIO.to_dict(SeqIO.parse(targets,"fasta"))
    illumina_adapters = SeqIO.to_dict(SeqIO.parse(adapters,"fasta"))
    snpdata = open(snpdata,"r").read().splitlines()[1:]
    snpdata = [x.split("\t") for x in snpdata]
    snpdata = {x[0]:[x[1],x[2],x[3]] for x in snpdata}


    pairs = {}
    for i in inread:
        if i.id not in pairs:
            pairs[i.id] = readpair()
            pairs[i.id].add1(name=i.id,r1=i)

    for i in inread2:
        if i.id in pairs:
            pairs[i.id].add2(name=i.id,r2=i)
        else:
            pairs[i.id] = readpair()
            pairs[i.id].add2(name=i.id,r2=i)

    for i in inoverlaps:
        pairs[i.id].add_ol(i)

    cutoff = 2

    for i in pairs:
        readpair_obj = pairs[i]

        readpair_obj.readcleanup(illumina_adapters)
        r1seq = readpair_obj.r1.seq
        r1seq = str(r1seq)
        r2seq = readpair_obj.r2.seq
        r2seq = str(r2seq)

        r1match = False
        r2match = False
        for primer in primers:
            p = primers[primer]
            primerseq = str(p.seq)
            rprimerseq = str(p.seq.reverse_complement())
            if not r1match:
                r1ham = hamming_distance(primerseq, r1seq[:len(primerseq) + 1])
                if r1ham <= cutoff:
                    r1match = p
            if not r1match:
                r1rcham = hamming_distance(rprimerseq, r1seq[:len(primerseq) + 1])
                if r1rcham <= cutoff:
                    r1match = SeqRecord.SeqRecord(p.seq.reverse_complement(),id="rc_"+p.id,name="",description="")
            if not r2match:
                r2ham = hamming_distance(primerseq, r2seq[:len(primerseq) + 1])
                if r2ham <= cutoff:
                    r2match = p
            if not r2match:
                r2rcham = hamming_distance(rprimerseq, r1seq[:len(primerseq) + 1])
                if r2rcham <= cutoff:
                    r2match = SeqRecord.SeqRecord(p.seq.reverse_complement(), id="rc_" + p.id, name="", description="")
        if not r1match:
            for primer in primers:
                p = primers[primer]
                primerseq = str(p.seq)
                r1swres = find_sw_match(primerseq, r1seq[:len(primerseq) + 1])
                if r1swres:
                    readstart = r1swres[0][0]
                    if readstart == 0:
                        r1match = p
        if not r2match:
            for primer in primers:
                p = primers[primer]
                primerseq = str(p.seq)
                r2swres = find_sw_match(primerseq, r2seq[:len(primerseq) + 1])
                if r2swres:
                    readstart = r2swres[0][0]
                    if readstart == 0:
                        r2match = p

        readpair_obj.add_primer_type(r1match,r2match)
        readpair_obj.check_products(nonprimertargets,targets,snpdata)
        readpair_obj.check_primer_dimer()

    outdetails = open(read_details,"w")
    outdetails.write("Name\tPairtype\tF_primer\tR_primer\tmixtype\tR1hasag\tR2hasag\tR1_prodmatch\tR2_prodmatch\toverlap_match\tR1_len\tR2_len\tollen\tprimerdimer\n") #self.name,self.type,self.p1,self.p2,self.mtype,self.r1hasag,self.r2hasag,self.r1len,self.r2len
    for i in pairs:
        readpair_obj = pairs[i]
        outdetails.write(readpair_obj.to_out_tab())
    outdetails.close()

    typecounts = {}
    correct_seqmatches = {}
    match_readlens = {}
    snpcalls = {}
    snpcounts = {}

    for i in pairs:
        readpair_obj = pairs[i]
        mtype = readpair_obj.mtype
        type = readpair_obj.type
        if type == "correct":
            match = False
            if readpair_obj.r1prodmatch and readpair_obj.r2prodmatch:
                match = True
            if readpair_obj.olmatch:
                match=True
            if mtype not in correct_seqmatches:
                correct_seqmatches[mtype] = 0
                if match:
                    correct_seqmatches[mtype] = 1
            else:
                if match:
                    correct_seqmatches[mtype] += 1
            if readpair_obj.snptype and match:
                snpinfo = readpair_obj.snptype
                gene,ref,snppos,snpref = snpinfo
                restuple = tuple(snpinfo)
                if gene not in snpcalls:
                    snpcalls[gene] = {restuple:1}
                else:
                    if restuple not in snpcalls[gene]:
                        snpcalls[gene][restuple] = 1
                    else:
                        snpcalls[gene][restuple]+=1



        if type not in typecounts:
            typecounts[type]={mtype:1}
        else:
            if mtype not in typecounts[type]:
                typecounts[type][mtype] = 1
            else:
                typecounts[type][mtype] += 1

        snpinfo = readpair_obj.snptype
        if snpinfo:
            if type not in snpcounts:
                snpcounts[type] = {mtype: {snpinfo:1}}
            else:
                if mtype not in snpcounts[type]:
                    snpcounts[type][mtype] = {snpinfo:1}
                else:
                    if snpinfo not in snpcounts[type][mtype]:
                        snpcounts[type][mtype][snpinfo] = 1
                    else:
                        snpcounts[type][mtype][snpinfo] += 1
        if type not in match_readlens:
            match_readlens[type] = {mtype:{"R1":[readpair_obj.r1len],"R2":[readpair_obj.r2len],"OL":[]}}
            if readpair_obj.olread:
                match_readlens[type][mtype]["OL"] = [len(readpair_obj.olread.seq)]

        else:
            if mtype not in match_readlens[type]:
                match_readlens[type][mtype] = {"R1":[readpair_obj.r1len],"R2":[readpair_obj.r2len],"OL":[]}
                if readpair_obj.olread:
                    match_readlens[type][mtype]["OL"] = [len(readpair_obj.olread.seq)]
            else:
                match_readlens[type][mtype]["R1"].append(readpair_obj.r1len)
                match_readlens[type][mtype]["R2"].append(readpair_obj.r2len)
                if readpair_obj.olread:
                    match_readlens[type][mtype]["OL"].append(len(readpair_obj.olread.seq))

    outf = open(output,"w")
    outf.write("Type\tR1\tR2\tcount\tinsert_correct\tR1_average_len\tR1_stdev\tR2_average_len\tR2_stdev\toverlap_avg\toverlap_stdev\tsnpcall\tsnptype\tsnperc\n")

    for type,mtyped in typecounts.items():
        for mtype,count in typecounts[type].items():
            if count > 0:
                r1avg = mean(match_readlens[type][mtype]["R1"])
                r2avg = mean(match_readlens[type][mtype]["R2"])
                if len(match_readlens[type][mtype]["OL"]) > 1:
                    olavg = mean(match_readlens[type][mtype]["OL"])
                    olstdev = stdev(match_readlens[type][mtype]["OL"])
                else:
                    olavg = 0
                    olstdev = 0
                if len(match_readlens[type][mtype]["R1"]) > 1:
                    r1stdev = stdev(match_readlens[type][mtype]["R1"])
                    r2stdev = stdev(match_readlens[type][mtype]["R2"])
                else:
                    r1stdev = 0
                    r2stdev = 0
                if mtype in correct_seqmatches:
                    matchno = correct_seqmatches[mtype]
                else:
                    matchno = 0
                if type in snpcounts:
                    if mtype in snpcounts[type]:
                        sndetails = snpcounts[type][mtype]
                        outcounts = {"ref": 0, "alt": 0, "other": 0}
                        nuc = {"ref": "", "alt": "", "other": ""}
                        for hit in sndetails:
                            outcounts[hit[1]] = sndetails[hit]
                            nuc[hit[1]] = hit[3]
                        majority = max(outcounts.values())
                        maxname = ""
                        for k, v in outcounts.items():
                            if v == majority:
                                maxname = k
                        majperc = (majority * 100) / (sum(outcounts.values()))
                        majnuc = nuc[maxname]
                    else:
                        majperc = ""
                        majnuc = ""
                        maxname = ""
                else:
                    majperc = ""
                    majnuc = ""
                    maxname = ""
                outf.write(f"{type}\t{mtype[0]}\t{mtype[1]}\t{count}\t{matchno}\t{r1avg:.2f}\t{r1stdev:.2f}\t{r2avg:.2f}\t{r2stdev:.2f}\t{olavg:.2f}\t{olstdev:.2f}\t{majnuc}\t{maxname}\t{majperc}\n")
    outf.close()


def main():

    r1 = sys.argv[1]
    r2 = sys.argv[2]
    olreads = sys.argv[3]
    primers = sys.argv[4]
    nonprimertargets = sys.argv[5]
    targets = sys.argv[6]
    adapters = sys.argv[7]
    output = sys.argv[8]
    snpdata = sys.argv[9]

    checkreads(r1,r2,olreads,primers,nonprimertargets,targets,adapters,output,snpdata)

if __name__ == '__main__':
    main()


