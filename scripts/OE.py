# 21-10-13
# kli
# mat -> OE/corr

import cooler
import numpy as np

def GetOE(mat):
    # Get Len 1.6m * index
    chr_len = mat.shape[0]
    # valid message < cutoff is gap!
    cut_off = chr_len/10
    # build mask
    mask = np.zeros(chr_len)
    num_mat = mat.copy()
    num_mat[num_mat > 0] = 1
    num_vector = np.sum(num_mat,axis=0)
    for i in range(chr_len):
        if(num_vector[i] >= cut_off):
            mask[i] = 1
    mask = mask == 1

    # get mean hash!
    ox = np.arange(chr_len)
    oy = np.arange(chr_len)
    omask = mask.copy()
    decay = {}
    for i in range(chr_len):
        o_diag = mat[(ox,oy)]
        o_diag_mask = o_diag[omask]
        # gap
        if(o_diag_mask.shape[0] == 0):
            decay[i] = 0
        else:
            decay[i] = o_diag_mask.mean()
        # left up.
        ox = np.delete(ox,-1)
        oy = np.delete(oy,0)
        omask = np.delete(omask,-1)

    # Get E mat!
    ex = np.arange(chr_len)
    ey = np.arange(chr_len)
    except_mat = np.ones_like(mat,dtype = np.float32)
    for i in range(chr_len):
        if(decay[i] == 0):
            # ex = np.delete(ex,-1)
            # ey = np.delete(ey,0)
            continue
        # xy <-> yx
        except_mat[(ex,ey)] = decay[i]
        except_mat[(ey,ex)] = decay[i]
        ex = np.delete(ex,-1)
        ey = np.delete(ey,0)
    
    # Get oe and corr
    oe = mat/except_mat
    cor_oe = np.corrcoef(np.cov(oe))
    cor_oe[np.isnan(cor_oe)] = 0

    return oe,cor_oe,mask


Pix1m = cooler.Cooler(snakemake.input[0]+"::resolutions/1600000")
mat = Pix1m.matrix(balance=False).fetch('chr1')
mat[np.isnan(mat)] = 0
oe,cor_oe,mask = GetOE(mat)
np.savetxt(snakemake.output[0],cor_oe,delimiter="\t")

