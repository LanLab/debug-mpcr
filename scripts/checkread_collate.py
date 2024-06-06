import sys


def checkread_collate(files,outprefix):
    files = list(files)
    pairtypes = []
    counts = {}
    intacts = {}
    strain_totals = {}
    inserts = {}
    snpcalls = {}
    for f in files:
        name = f.replace("_paircounts.txt","")
        counts[name] = {}
        intacts[name] = {}
        inserts[name] = {}
        snpcalls[name] = {}
        strain_totals[name] = 0
        inp = open(f,"r").read().splitlines()
        for line in inp[1:]:
            col = line.split("\t")
            pair = (col[0],col[1],col[2])
            if pair not in pairtypes:
                pairtypes.append(pair)
            counts[name][pair] = col[3]
            intacts[name][pair] = col[4]
            inserts[name][pair] = col[9]
            strain_totals[name] += int(col[3])
            snpcalls[name][pair] = (col[11],col[12],col[13])
    out = open(f"{outprefix}_primer_pair_results.txt","w")
    outperc = open(f"{outprefix}_primer_pair_results_perc.txt","w")
    out2 = open(f"{outprefix}_primer_pair_results_w_intact_insert.txt","w")
    out2perc = open(f"{outprefix}_primer_pair_results_w_intact_insert_perc.txt","w")
    out3 =open(f"{outprefix}_primer_pair_results_inserts.txt","w")
    sampleorder = list(sorted(counts.keys()))#,key=lambda x:x.split("_")[1]))
    out.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    outperc.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    outsnp = open(f"{outprefix}_primer_pair_snp_typing.txt","w")
    outsnp2 = open(f"{outprefix}_primer_pair_snp_typing_details.txt", "w")
    outperc.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    for pair in pairtypes:
        pairls = list(pair)
        pairls.append(":".join(list(pair)))
        pairperc1ls = list(pair)
        pairperc1ls.append(":".join(list(pair)))
        for strain in sampleorder:
            if pair in counts[strain]:
                count = counts[strain][pair]
                pairls.append(count)
                pairperc1ls.append("{:.4f}".format(float(100 * int(count)) / int(strain_totals[strain])))
            else:
                pairls.append('0')
                pairperc1ls.append('0')
        out.write("\t".join(pairls)+"\n")
        outperc.write("\t".join(pairperc1ls) + "\n")
    out.close()
    outperc.close()

    out2.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    out2perc.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    for pair in pairtypes:
        pairls = list(pair)
        pairpercls = list(pair)
        pairpercls.append(":".join(list(pair)))
        for strain in sampleorder:
            if pair in intacts[strain]:
                pairls.append(intacts[strain][pair])
                pairpercls.append("{:.2f}".format(float(100*int(intacts[strain][pair]))/int(counts[strain][pair])))
            else:
                pairls.append('0')
                pairpercls.append('0')
        out2.write("\t".join(pairls)+"\n")
        out2perc.write("\t".join(pairpercls)+"\n")
    out2.close()
    out2perc.close()

    out3.write("pairtype\tread1_primer\tread2_primer\tmixtype\t{}\n".format("\t".join(sampleorder)))
    for pair in pairtypes:
        pairls = list(pair)
        pairls.append(":".join(list(pair)))
        for strain in sampleorder:
            if pair in inserts[strain]:
                count = inserts[strain][pair]
                pairls.append(count)
            else:
                pairls.append('0')
        out3.write("\t".join(pairls)+"\n")

    out3.close()

    outsnp.write("gene\t{}\n".format("\t".join(sampleorder)))
    outsnp2.write("gene\t{}\n".format("\t".join(sampleorder)))
    for pair in pairtypes:
        if pair[0] == "correct":
            gene = pair[1].replace("_f","")
            pairls = [gene]
            pairls2 = [gene]
            for strain in sampleorder:
                if pair in snpcalls[strain]:
                    nuccall = snpcalls[strain][pair]
                    pairls.append(nuccall[1])
                    pairls2.append(",".join(list(nuccall)))
                else:
                    pairls.append('no_reads')
                    pairls2.append('no_reads')
            outsnp.write("\t".join(pairls)+"\n")
            outsnp2.write("\t".join(pairls2)+"\n")

    outsnp.close()

def main():
    files = sys.argv[2:]
    outprefix = sys.argv[1]
    checkread_collate(files,outprefix)

if __name__ == '__main__':
    main()