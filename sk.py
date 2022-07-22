# 10-8
# kli
# zero flow
# ...

import re
import os

# indir:in  seqformate:{spec}_R*.fq
def GetKeys():
    seq_path = os.getcwd() + '/' + 'in'
    seq_names = os.listdir(seq_path)
    key_words = []
    for sn in seq_names:
        ref = re.match("(.*?)_.*fq",sn)
        key_words.append(ref.groups()[0])
    return  list(set(key_words))

keys = GetKeys()

rule all:
    input:
        expand("cool/{spc}.mcool",spc = keys),
        expand("about/{spc}.txt",spc = keys)
        # expand("tadout/{spc}_rect.pdf",spc = keys),
        # "out/go1.pdf",
        # "out/go2.pdf"

rule alig:
    input:
        r1 = "in/{spec}_R1.fq",
        r2 = "in/{spec}_R2.fq",
        rf = "refer/ref.fa"
    output:
        "mapped/{spec}.bam"
    shell:
        "bwa mem -SP5M -t8 {input.rf} {input.r1} {input.r2} | samtools "
        "view -bhS - > {output}"

rule Parse:
    input:
        b = "mapped/{pin}.bam",
        ch = "refer/chrom.c"
    output:
        "mapped/{pin}.pairsam"
    shell:
        "samtools view -h {input.b} | pairtools parse -c {input.ch} -o {output}"

rule Sort:
    input:
        "mapped/{sin}.pairsam"
    output:
        "mapped/{sin}.sorted.pairsam"
    shell:
        "pairtools sort --nproc 8 -o {output} {input}"

rule Dup:
    input:
        "mapped/{din}.sorted.pairsam"
    output:
        "mapped/{din}.duped.pairsam"
    shell:
        "pairtools dedup --mark-dups -o {output} {input}"

rule Select:
    input:
        "mapped/{selin}.duped.pairsam"
    output:
        "mapped/{selin}.filtered.pairsam"
    shell:
        "pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\")"
        " or (pair_type == \"RU\")' -o {output} {input}"

rule Split:
    input:
        "mapped/{spin}.filtered.pairsam"
    output:
        "cool/{spin}.pairs"
    shell:
        "pairtools split --output-pairs {output} {input}"

# !!! gunzip <-> bgzip
# under index we could...
# pairix out.pairs.gz "chr1:1-60000|chr4:300-50000"
rule Index:
    input:
        "cool/{inin}.pairs"
    output:
        "cool/{inin}.pairs.gz"
    shell:
        "bgzip {input} && pairix -f {output}"


# cooler dump -t pixels --header --join -r chr19 out.cool
rule Getcool:
    input:
        gp = "cool/{gin}.pairs.gz",
        ch = "refer/chrom.c"
    output:
        "cool/{gin}.cool"
    shell:
        "cooler cload pairix {input.ch}:50000 {input.gp} {output}"

# Norm is in zoom step.
# after zoom,we dump will use "out.mcool::resolution/50000"
# zoom is x 2x 4x .. 2^nx
rule Zoom:
    input:
        "cool/{zin}.cool"
    output:
        "cool/{zin}.mcool"
    shell:
        "cooler balance {input} && cooler zoomify --balance -o "
        "{output} {input}"


# Tad Part
# From there ,we input .mcool ,
# but in our .py we generate .mat and then 
# we use OnTAD to generate a .mat.tad file.
# bue snk can't identity it,so we pretend generate .mat!
# we cheat the snakemake!
# There will generate the .tad file.
# rule GetMarix:
    # input:
        # "cool/{gin}.mcool"
    # output:
        # "tad/{gin}.mat"
    # script:
        # "scripts/GetMar.py"
# pass flag! Real tad genetare is in GetMarix!!!
# rule GetTad:
    # input:
        # "tad/{tin}.mat"
    # output:
        # "tad/{tin}.tad"
    # script:
        # "scripts/Tad.py"
# rule TadPlot:
    # input:
        # 'tad/{tpin}.tad'
    # output:
        # "tadout/{tpin}_line.pdf",
        # "tadout/{tpin}_rect.pdf"
    # script:
        # "scripts/TadPlot.R"


# AB part
rule GetOE:
    input:
        "cool/{oein}.mcool"
    output:
        "about/{oein}.txt"
    script:
        "scripts/OE.py"
# rule ABPlot:
    # input:
        # "about/{abin}.txt"
    # output:
        # "about/{abin}.pdf"
    # script:
        # "scripts/CompartmentPlot.R"

# Loop Part

# Group Part
# rule gr1:
    # input:
        # expand("cool/{ain}.mcool",ain = keys)
    # output:
        # "out/go1.pdf"
    # script:
        # "scripts/group1.py"

# rule gt2:
    # input:
        # expand("cool/{ain}.mcool",ain = keys)
    # output:
        # "out/go2.pdf"
    # script:
        # "scripts/group2.py"