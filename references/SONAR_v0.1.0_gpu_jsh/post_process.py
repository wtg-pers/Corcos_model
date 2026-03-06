# import for system.
import os
import sys
from sys import platform
# import for numerical calculation.
import numpy as np
import cupy as cp
import math
from numba import jit
from itertools import product
from scipy.fftpack import fft
from scipy.signal import welch
from scipy.signal.windows import hann

# lbn_ob3 = [17.8, 22.4, 28.2, 35.5, 44.7, 56.2, 70.8, 89.1, 112.0,
#            141.0, 178.0, 224.0, 282.0, 355.0, 447.0, 562.0, 708.0,
#            891.0, 1120.0, 1410.0, 1780.0, 2240.0, 2820.0, 3550.0,
#            4470.0, 5620.0, 7080.0, 8910.0, 11200.0, 14130.0, 17780.0]
# cen_ob3 = [20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0,
#            160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0,
#            1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0,
#            5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0,
#            20000.0]
# ubn_ob3 = [22.4, 28.2, 35.5, 44.7, 56.2, 70.8, 89.1, 112.0, 141.0,
#            178.0, 224.0, 282.0, 355.0, 447.0, 562.0, 708.0, 891.0,
#            1120.0, 1410.0, 1780.0, 2240.0, 2820.0, 3550.0, 4470.0,
#            5620.0, 7080.0, 8910.0, 11200.0, 14130.0, 17780.0, 22390.0]

# fl = lbn_ob3
# fc = cen_ob3
# fu = ubn_ob3
# lfrq = lbn_ob3
# cfrq = cen_ob3
# ufrq = ubn_ob3


def delete_previous_FFT(obj_par):
    NBLADE, RPM = obj_par.NBLADE, int(obj_par.RPM)
    fft_out = obj_par.fft_out
    # Set deletion command for the current system platform.
    if platform == 'linux' or platform == 'linux2' or platform == 'darwin':
        def delete_file(file_path):
            os.system(f'rm -f {file_path}')
    elif platform == 'win32':
        def delete_file(file_path):
            os.system(f'del /q {file_path}')

    # Run the deletion operation corresponding to the output directory,
    if os.path.isdir(fft_out):  # when the directory exists.
        ton_check = '%sb_%sRPM_tonal.dat' % (NBLADE, RPM)
        for filename in os.listdir(fft_out):
            if ton_check == filename:
                file_path = os.path.join(fft_out, filename)
                delete_file(file_path)
                print(f"{filename} is deleted.")

        bbn_check = '%sb_%sRPM_bbn.dat' % (NBLADE, RPM)
        for filename in os.listdir(fft_out):
            if bbn_check == filename:
                file_path = os.path.join(fft_out, filename)
                delete_file(file_path)
                print(f"{filename} is deleted.")

        oaspl_check = '%sb_%sRPM_oaspl.dat' % (NBLADE, RPM)
        for filename in os.listdir(fft_out):
            if oaspl_check == filename:
                file_path = os.path.join(fft_out, filename)
                delete_file(file_path)
                print(f"{filename} is deleted.")
    else:  # when the directory is newly created.
        os.system('mkdir %s' % fft_out)
    return None


def a_weighting(f):
    c1 = 12194**2
    c2 = 20.6**2
    c3 = 107.7**2
    c4 = 737.9**2

    A1000 = ((c1 * 1000**4) /
             ((1000**2 + c2) * np.sqrt((1000**2 + c3) *
                                       (1000**2 + c4)) * (1000**2 + c1)))

    num = c1 * f**4
    den = (f**2 + c2) * np.sqrt((f**2 + c3) * (f**2 + c4)) * (f**2 + c1)

    ra = num / den

    a = 20 * np.log10(ra) - 20 * np.log10(A1000)

    return a


# @jit(nopython=True, cache=True)
def convert_nb_to_ob_data(lfrq, cfrq, ufrq, nb_frq, nb_spl):
    # Get necessary variables.
    nfn = len(nb_frq)  # Number of NB frequencies
    nfo = len(cfrq)  # Number of OB center frequencies

    # Initialize variables.
    dat_ob = np.zeros(nfo)
    dat_oaspl = 0.0

    # Compute df (this is a common value for the overall frequency range).
    df = nb_frq[1] - nb_frq[0]

    for j in range(nfo):
        for i in range(nfn):
            frac = 0.0
            nb_lfrq = nb_frq[i] - 0.5 * df
            nb_ufrq = nb_frq[i] + 0.5 * df

            # Filter unnecessary computation.
            if (nb_ufrq < lfrq[j]) or (ufrq[j] < nb_lfrq):
                continue

            if (lfrq[j] < nb_lfrq) and (nb_ufrq < ufrq[j]):
                frac = nb_ufrq - nb_lfrq
            elif (lfrq[j] < nb_lfrq) and (nb_ufrq >= ufrq[j]):
                frac = ufrq[j] - nb_lfrq
            elif (lfrq[j] >= nb_lfrq) and (nb_ufrq < ufrq[j]):
                frac = nb_ufrq - lfrq[j]
            elif (lfrq[j] >= nb_lfrq) and (nb_ufrq >= ufrq[j]):
                frac = ufrq[j] - lfrq[j]

            dat_ob[j] += (1.0/df) * frac * 10.0**(nb_spl[i]/10.0)

    # Set zero SPL values,
    for f in range(0, nfo):
        if np.isnan(dat_ob[f]):  # when the value is nan.
            dat_ob[f] = 1.0
            continue
        if np.isinf(dat_ob[f]):  # when the value is infinity.
            dat_ob[f] = 1.0
            continue
        if (dat_ob[f] <= 0.0):  # when the value is negative.
            dat_ob[f] = 1.0
            continue
        dat_oaspl += dat_ob[f]

    dat_ob = 10.0 * np.log10(dat_ob)
    if (dat_oaspl >= 1.0):
        dat_oaspl = 10.0 * np.log10(dat_oaspl)
    else:
        dat_oaspl = 0.0

    return dat_ob, dat_oaspl


def do_FFT(time_data, p_data, df):
    N = len(time_data)
    fs = N / time_data[-1]

    seg_leng = int(fs // df)
    noverlap = seg_leng // 2  # 50% overlap

    step = seg_leng - noverlap
    seg_idx = range(0, len(p_data) - seg_leng + 1, step)

    fft_results = []
    for start in seg_idx:
        segment = p_data[start:start+seg_leng]
        nnyq = (len(segment) / 2)
        window = np.hanning(len(segment))
        beta = np.sum(window**2) / len(segment)
        sig = np.multiply(segment, window)
        fft_rst = fft(sig)
        fft_amp = np.abs(fft_rst)
        fft_ang = np.angle(fft_rst, deg=False)
        fft_psd = fft_amp**2.0 / (beta * nnyq)
        amp = np.sqrt(fft_psd / nnyq) / np.sqrt(2)
        amp_trim = amp[:int(nnyq)]
        ang_trim = fft_ang[:int(nnyq)]

        nb_frq = (fs / 2) * np.linspace(0, 1, int(nnyq))

        fft_results.append(amp_trim)

    fft_avg = np.mean(np.abs(fft_results), axis=0)
    nb_spl = 20 * np.log10(fft_avg / 2e-05)

    return nb_frq, nb_spl


def do_FFT_processing(obj_par, obj_rot, obj_mic, obj_sig):
    NBLADE, RPM = obj_par.NBLADE, int(obj_par.RPM)
    nmics = obj_par.nmics
    mics = cp.asnumpy(obj_mic.mics)
    map_x = mics[:, 0]
    map_y = mics[:, 1]
    map_z = mics[:, 2]

    hxy = []
    r = []
    el = []
    az = []
    for i in range(len(map_x)):
        hxy_temp = math.hypot(map_x[i], map_y[i])
        r_temp = math.hypot(hxy_temp, map_z[i])
        el_temp = math.atan2(map_z[i], hxy_temp)
        az_temp = math.atan2(map_y[i], map_x[i])

        hxy.append(hxy_temp)
        r.append(r_temp)
        el.append(el_temp)
        az.append(az_temp)

    el = np.rad2deg(el)

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        lfrq = cp.asnumpy(obj_sig.lbn_ob)
        cfrq = cp.asnumpy(obj_sig.cen_ob)
        ufrq = cp.asnumpy(obj_sig.ubn_ob)
    elif FRQ_BAND == 'third':
        lfrq = cp.asnumpy(obj_sig.lbn_ob3)
        cfrq = cp.asnumpy(obj_sig.cen_ob3)
        ufrq = cp.asnumpy(obj_sig.ubn_ob3)
    elif FRQ_BAND == 'twelfth':
        lfrq = cp.asnumpy(obj_sig.lbn_ob12)
        cfrq = cp.asnumpy(obj_sig.cen_ob12)
        ufrq = cp.asnumpy(obj_sig.ubn_ob12)

    fft_out = obj_par.fft_out

    cp_new_t2 = obj_rot.new_t2
    cp_ld_p2, cp_tk_p2 = obj_rot.ld_p2, obj_rot.tk_p2
    cp_tonal_p = obj_rot.tonal_p

    p8 = obj_rot.p8  # Airfoil-Self noise
    p9 = obj_rot.p9  # BWI noise
    p10 = obj_rot.p10  # Airfoil-Self noise + BWI noise

    new_t2 = cp.asnumpy(cp_new_t2)
    ld_p2 = cp.asnumpy(cp_ld_p2)
    tk_p2 = cp.asnumpy(cp_tk_p2)
    tonal_p = cp.asnumpy(cp_tonal_p)

    p8 = cp.asnumpy(p8)
    p9 = cp.asnumpy(p9)
    p10 = cp.asnumpy(p10)

    nb_frq_ld = []
    nb_frq_tk = []
    nb_frq = []

    nb_spl_ld = []
    nb_spl_tk = []
    nb_spl = []

    ovlp = 0.5

    for i in range(nmics):
        tmp_idx = []
        tmp_ob = []
        tmp_nb_ld = []
        tmp_nb_tk = []
        tmp_nb = []

        out_tonal = os.path.join(fft_out, '%sb_%sRPM_tonal.dat' % (
            NBLADE, RPM))
        out_bbn = os.path.join(fft_out, '%sb_%sRPM_bbn.dat' % (
            NBLADE, RPM))
        out_oaspl = os.path.join(fft_out, '%sb_%sRPM_oaspl.dat' % (
            NBLADE, RPM))

        time_data = new_t2[i, :]
        ld_data = ld_p2[i, :]
        tk_data = tk_p2[i, :]
        tonal_data = tonal_p[i, :]

        df = 3.125
        # df = 5
        # df = 1

        nb_frq_ld, nb_spl_ld = do_FFT(time_data, ld_data, df)
        nb_frq_tk, nb_spl_tk = do_FFT(time_data, tk_data, df)
        nb_frq, nb_spl = do_FFT(time_data, tonal_data, df)

        ob_spl, ob_oaspl = convert_nb_to_ob_data(lfrq, cfrq, ufrq,
                                                 nb_frq, nb_spl)

        spl_bbn_as = 10 * np.log10(10**(p8[i, :] / 10))
        spl_bbn_bwi = 10 * np.log10(10**(p9[i, :] / 10))
        spl_bbn_total = 10 * np.log10(10**(ob_spl/10) + 10**(p10[i, :] / 10))
        spl_tonal_as = 10 * np.log10(10**(ob_spl/10) + 10**(p8[i, :] / 10))
        spl_tonal_bwi = 10 * np.log10(10**(ob_spl/10) + 10**(p9[i, :] / 10))

        oaspl_as = 10 * np.log10(np.sum(10**(spl_bbn_as / 10)))
        oaspl_bwi = 10 * np.log10(np.sum(10**(spl_bbn_bwi / 10)))
        oaspl_total = 10 * np.log10(np.sum(10**(spl_bbn_total / 10)))

        a_correction = a_weighting(cfrq)
        a_weighted_spl = spl_bbn_total + a_correction
        # a_linear_spl = np.sum(10 ** (a_weighted_spl / 10))
        # spl_a = 10 * np.log10(a_linear_spl)
        spl_a = a_weighted_spl

        f1 = open(out_tonal, 'a')
        f2 = open(out_bbn, 'a')
        f3 = open(out_oaspl, 'a')
        if i == 0:
            f1.write('variables="frequency[Hz]", '
                     '"SPL<sub>Loading</sub>[dB]", '
                     '"SPL<sub>Thickness</sub>[dB]", '
                     '"SPL<sub>Loading+Thickness</sub>[dB]"\n')

            f2.write('variables="1/3 Octave band frequency[Hz]", '
                     '"SPL<sub>tonal</sub>[dB]", '
                     '"SPL<sub>AS</sub>[dB]", '
                     '"SPL<sub>bwi</sub>[dB]", '
                     '"SPL<sub>tonal+AS</sub>[dB]", '
                     '"SPL<sub>tonal+BWI</sub>[dB]", '
                     '"SPL<sub>tonal+AS+BWI</sub>[dB]", '
                     '"SPL<sub>total</sub>[dB(A)]"\n')

            f3.write('variables="Directivity[deg]", '
                     '"OASPL[dB]<sub>tonal</sub>", '
                     '"OASPL[dB]<sub>AS</sub>", '
                     '"OASPL[dB]<sub>BWI</sub>", '
                     '"OASPL[dB]<sub>tonal+AS+BWI</sub>"\n')
            f3.write(f'zone T="pressent_{NBLADE}b_{RPM}RPM"\n')

        # elif i == (nmics-1):
        #     f3.write('zone T="near field"\n')

        f1.write(f'zone t= "Mic{i+1}"\n')

        f2.write(f'zone t= "Mic{i+1}"\n')

        for ii in range(len(nb_spl)):
            f1.write(f'{nb_frq[ii]}\t'
                     f'{nb_spl_ld[ii]}\t'
                     f'{nb_spl_tk[ii]}\t'
                     f'{nb_spl[ii]}\n')

        for jj in range(len(cfrq)):
            f2.write(f'{cfrq[jj]}\t'
                     f'{ob_spl[jj]}\t'
                     f'{p8[i, jj]}\t'
                     f'{p9[i, jj]}\t'
                     f'{spl_tonal_as[jj]}\t'
                     f'{spl_tonal_bwi[jj]}\t'
                     f'{spl_bbn_total[jj]}\t'
                     f'{spl_a[jj]}\n')

        f3.write(f"{el[i]:.2f}\t"
                 f"{ob_oaspl:.2f}\t"
                 f"{oaspl_as:.2f}\t"
                 f"{oaspl_bwi:.2f}\t"
                 f"{oaspl_total:.2f}\n")

    f1.close()
    f2.close()
    f3.close()

    return None


def do_process(objs):
    # Get variables.
    obj_par = objs[0]
    obj_rot = objs[1]
    obj_mic = objs[2]
    obj_sig = objs[3]

    delete_previous_FFT(obj_par)
    do_FFT_processing(obj_par, obj_rot, obj_mic, obj_sig)

    return None
