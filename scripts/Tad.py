
with open(snakemake.input[0]+".tad","r") as f:
    a = f.read()

with open(snakemake.output[0],"w") as f:
    f.write(a)
