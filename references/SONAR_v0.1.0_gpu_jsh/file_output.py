# import for system.
import os
# import for numerical calculation.
import numpy as np
import cupy as cp
from tqdm import tqdm


def write_tonal_noise(obj_par, obj_rot, obj_mic):
    NBLADE, RPM = obj_par.NBLADE, int(obj_par.RPM)
    current_file = os.path.join(obj_par.dir_out, "%sb_%srpm_tonal_p.dat" % (
        NBLADE, RPM))
    nmics = obj_par.nmics
    NINTERP = obj_par.NINTERP
    period = obj_par.dt * float(obj_par.NSTEP)

    cp_mics = obj_mic.mics

    cp_new_t2 = obj_rot.new_t2
    cp_ld_p2, cp_tk_p2 = obj_rot.ld_p2, obj_rot.tk_p2
    cp_tonal_p = obj_rot.tonal_p

    mics = cp.asnumpy(cp_mics)

    new_t2 = cp.asnumpy(cp_new_t2)
    ld_p2 = cp.asnumpy(cp_ld_p2)
    tk_p2 = cp.asnumpy(cp_tk_p2)
    tonal_p = cp.asnumpy(cp_tonal_p)

    with open(current_file, 'w') as f:
        for idx in range(nmics):
            f.write('variables="time","revolution",'
                    '"loading pressure","thickness pressure",'
                    '"tonal pressure"\n')
            f.write(f'zone T="mic#{idx + 1}, '
                    f'coord={mics[idx, 0]}, {mics[idx, 1]}, '
                    f'{mics[idx, 2]}"\n')
            for j in range(NINTERP):
                tt = new_t2[idx, j]
                f.write(f"{tt}\t{tt/period}\t{ld_p2[idx, j]}\t"
                        f"{tk_p2[idx, j]}\t{tonal_p[idx, j]}\n")

    return None


def write_broadband_noise(obj_par, obj_rot, obj_mic, obj_sig):
    NBLADE, RPM = obj_par.NBLADE, int(obj_par.RPM)
    current_file = os.path.join(obj_par.dir_out, "%sb_%srpm_broadband_spl.dat"
                                % (NBLADE, RPM))
    nmics = obj_par.nmics
    nfreq = obj_par.nfreq

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    cp_mics = obj_mic.mics
    mics = cp.asnumpy(cp_mics)

    p1, p2, p3, p4 = obj_rot.p1, obj_rot.p2, obj_rot.p3, obj_rot.p4
    p5, p6, p7, p8 = obj_rot.p5, obj_rot.p6, obj_rot.p7, obj_rot.p8
    p9 = obj_rot.p9
    p10 = obj_rot.p10
    p1, p2, p3 = cp.asnumpy(p1), cp.asnumpy(p2), cp.asnumpy(p3)
    p4, p5, p6 = cp.asnumpy(p4), cp.asnumpy(p5), cp.asnumpy(p6)
    p7, p8 = cp.asnumpy(p7), cp.asnumpy(p8)
    p9 = cp.asnumpy(p9)
    p10 = cp.asnumpy(p10)

    with open(current_file, 'w') as f1:
        for i in range(nmics):
            f1.write('variables= "frequency[Hz]", "TBL_TE<sub>p</sub>[dB]", '
                     '"TBL_TE<sub>s</sub>[dB]", "TBL_TE<sub>AoA</sub>[dB]", '
                     '"TBL_TE<sub>total</sub>[dB]", "LBL<sub>VS</sub>[dB]", '
                     '"TE_BLUNTNESS[dB]", "TIP<sub>VS</sub>[dB]", '
                     '"Airfoil_Self_NOISE<sub>total</sub>[dB]", "BWI[dB]", '
                     '"Broadband Noise[dB]"\n')
            f1.write(f'zone T="mic#{i + 1}, '
                     f'coord={mics[i, 0]}, {mics[i, 1]}, '
                     f'{mics[i, 2]}"\n')
            for j in range(nfreq):
                f1.write(f"{cen_ob[j]}\t"
                         f"{p1[i, j]}\t"
                         f"{p2[i, j]}\t"
                         f"{p3[i, j]}\t"
                         f"{p4[i, j]}\t"
                         f"{p5[i, j]}\t"
                         f"{p6[i, j]}\t"
                         f"{p7[i, j]}\t"
                         f"{p8[i, j]}\t"
                         f"{p9[i, j]}\t"
                         f"{p10[i, j]}\n"
                         )

    return None
