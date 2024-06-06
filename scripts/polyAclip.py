from Bio import SeqIO
from Bio import SeqIO,SeqRecord,Seq
import sys
from time import sleep as sl


def remove_polya(r):
    polya = "AAAA"
    polyg = "GGGGGGG"
    polyt = "TTTTTTT"
    rseq = str(r.seq)

    apos = rseq.rfind(polya)
    gpos = rseq.find(polyg)
    tpos = rseq.find(polyt)
    if apos > 0 and (gpos > 0 or tpos > 0):

        newrseq = rseq[:apos]
        r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                      name=r.name, letter_annotations={
                "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})

    elif gpos > 0 and "GCCGCTGGCGCCCGCGATGGCGATGC" not in rseq:
        newrseq = rseq[:gpos]
        r = SeqRecord.SeqRecord(Seq.Seq(newrseq), id=r.id, description=r.description,
                                name=r.name, letter_annotations={
                "phred_quality": r.letter_annotations["phred_quality"][:len(newrseq)]})


    return r

def main():
    inf = sys.argv[1]
    outf = sys.argv[2]
    inreads = SeqIO.parse(inf,"fastq")
    outreads = []
    c = 0
    for read in inreads:
        trimread = remove_polya(read)
        if len(trimread.seq) != len(read.seq):
            c+=1
        outreads.append(trimread)
    print(c)
    SeqIO.write(outreads,outf,"fastq")

if __name__ == '__main__':
    main()