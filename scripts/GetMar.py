import cooler
import numpy as np
import os

Pix50k = cooler.Cooler(snakemake.input[0]+"::resolutions/50000")
mat = Pix50k.matrix(balance=True).fetch('chr21')
mat[np.isnan(mat)] = 0
mat = np.around(mat,decimals=3)

np.savetxt(snakemake.output[0],mat,delimiter="\t",fmt="%.3f")

# This file could change the args.
os.system("OnTAD "+snakemake.output[0]+" -o "+snakemake.output[0])
