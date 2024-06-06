

configfile: "config.json"
pref=config["inprefix"]
suff=config["insuffix"]
outpref = config["outpref"]

CONDITIONS = glob_wildcards(pref+"{sample}"+suff).sample

conda environment.yml

rule all:
    input: outpref+"_primer_pair_results.txt"
    default_target: True

checkpoint gunzip:
    input: pref+"{sample}_R{fr}_001.fastq.gz"
    output: "{sample}_R{fr}_001.fastq"
    group:"group1"
    conda:
        "mpcr_testing"
    shell: """
    gunzip -c {input} > {output}
    """

checkpoint runclip:
    input: "{sample}_R{fr}_001.fastq"
    output: "{sample}_R{fr}_001_clip.fastq"
    group:"group1"
    conda:
        "mpcr_testing"
    shell:"""
    python scripts/polyAclip.py {input} {output}
    """

checkpoint runfastp:
    input:
        in1="{sample}_R1_001_clip.fastq",
        in2="{sample}_R2_001_clip.fastq"
    output:
        out1="{sample}_R1_001_clip_fastp.fastq",
        out2="{sample}_R2_001_clip_fastp.fastq",
    group:"group2"
    conda:
        "mpcr_testing"
    shell:"""
    fastp -c -i {input.in1} -o {output.out1} -I {input.in2} -O {output.out2}
    rm fastp.html
    rm fastp.json
    """

checkpoint runflash:
    input:
        in1="{sample}_R1_001_clip_fastp.fastq",
        in2="{sample}_R2_001_clip_fastp.fastq"
    output:
        output="{sample}_OL_001_clip_fastp_flash_extendedFrags.fastq"
    group:"group2"
    conda:
        "mpcr_testing"
    params:
       pfx ="{sample}_OL_001_clip_fastp_flash"
    shell:"""
    flash --max-overlap 150 --allow-outies {input.in1} {input.in2} -o {params.pfx}
    mv {params.pfx}.extendedFrags.fastq {params.pfx}_extendedFrags.fastq
    rm {params.pfx}.*
    """

checkpoint runcheckreads:
    input:
        in1 = "{sample}_R1_001_clip_fastp.fastq",
        in2 = "{sample}_R2_001_clip_fastp.fastq",
        inol="{sample}_OL_001_clip_fastp_flash_extendedFrags.fastq"
    output: "{sample}_paircounts.txt"
    group:"group2"
    conda:
        "mpcr_testing"
    params:
       primers = config["primers"],
       nonprimertargets = config["nonprimertargets"],
       targets = config["targets"],
       adapters = config["adapters"],
       snpdata = config["snpdata"]
    shell:"""
        python scripts/checkreads.py {input.in1} {input.in2} {input.inol} {params.primers} {params.nonprimertargets} {params.targets} {params.adapters} {output} {params.snpdata}
    """




rule collated_summary:
    input: expand("{group}_paircounts.txt",group=CONDITIONS)
    output: outpref+"_primer_pair_results.txt"
    group: "group3"
    params:
       pfx =outpref
    conda:
        "mpcr_testing"
    shell:"""
    python scripts/checkread_collate.py {params.pfx} {input}
    """
