# import for system.
import os
import time
# import for numerical calculation.
import math
import numpy as np
import cupy as cp
from tqdm import tqdm

# import for methods and functions.
import file_input as inp
import geometry as geo
import noise_calculation as calc
import noise_model as model
import file_output as out


def do_process(i1, args, logs, objs):
    """

    Parameters
    ----------
    i1 : int
        First index of process.
    args : class
        Arguments object.
    logs : class
        loggers object.
    objs : list
        list of objects.

    Returns
    -------
    None.

    """

    # Get variables.
    obj_par = objs[0]
    obj_rot = objs[1]
    obj_mic = objs[2]
    obj_sig = objs[3]
    MODE = obj_par.MODE
    TON_OPTIONS = obj_par.TON_OPTIONS
    FRQ_BAND = obj_par.FRQ_BAND

    # toggle options.
    NOISE_TON = obj_par.NOISE_TON
    NOISE_BBN = obj_par.NOISE_BBN
    NOISE_LAM = obj_par.NOISE_LAM
    NOISE_TUR = obj_par.NOISE_TUR
    NOISE_BLT = obj_par.NOISE_BLT
    NOISE_TIP = obj_par.NOISE_TIP
    NOISE_BWI = obj_par.NOISE_BWI

    nmics = obj_par.nmics
    OM, ON = obj_par.OM, obj_par.ON
    NINTERP = obj_par.NINTERP
    NSTEP, NBLADE, NSEG = obj_par.NSTEP, obj_par.NBLADE, obj_par.NSEG

    temp_mics = obj_mic.temp_mics

    new_t2 = obj_rot.new_t2
    ld_p2 = obj_rot.ld_p2
    tk_p2 = obj_rot.tk_p2
    tonal_p = obj_rot.tonal_p

    p1, p2, p3, p4 = obj_rot.p1, obj_rot.p2, obj_rot.p3, obj_rot.p4
    p5, p6, p7, p8 = obj_rot.p5, obj_rot.p6, obj_rot.p7, obj_rot.p8
    p9, p10 = obj_rot.p9, obj_rot.p10

    # Initialize and set variables.
    i2 = 0

    test_f = open('lti_temp.dat', 'w')
    test_f.write('variables="Updated Index of k<sup>*</sup><sub>2</sub>", '
                 '"Microphone Number", '
                 '"Magnitude of L<sub>TI</sub>"\n')
    test_f.write(f'zone t="1", i=2, j=17, k=14\n')
    hxy = []
    r = []
    el = []
    az = []

    if MODE != 'post':
        i3 = 0
        # i2 += 1
        # logs.info(' --- %i.%i    Reading microphone coordinate datas.' %
        #           (i1, i2))
        # inp.get_microphone_coordinates(obj_par, obj_mic)
        mics = obj_mic.mics
        obs_x = mics[:, 0]
        obs_y = mics[:, 1]
        obs_z = mics[:, 2]

        for i in range(len(obs_x)):
            hxy_temp = math.hypot(obs_x[i], obs_y[i])
            r_temp = math.hypot(hxy_temp, obs_z[i])
            el_temp = math.atan2(obs_z[i], hxy_temp)
            az_temp = math.atan2(obs_y[i], obs_x[i])

            hxy.append(hxy_temp)
            r.append(r_temp)
            el.append(el_temp)
            az.append(az_temp)

        el = np.rad2deg(el)

        # Aerodynamic reference value.
        mach = (obj_rot.ufree / obj_par.C0).astype(cp.float64)
        re_chord = (obj_rot.ufree * obj_rot.C * obj_par.RHO /
                    obj_par.MU).astype(cp.float64)
        obj_rot.re_chord, obj_rot.mach = re_chord, mach

        i2 += 1
        logs.info(' --- %i.%i    Propagating surface data to microphone.' %
                  (i1, i2))
        # ts_noise = time.time()
        for m in tqdm(range(nmics), desc="    "):
            temp_mics[:] = cp.array([obs_x[m], obs_y[m], obs_z[m]])
            # Initialize the pressure component of broadband.
            p1[m, :] = 0.0
            p2[m, :] = 0.0
            p3[m, :] = 0.0
            p4[m, :] = 0.0
            p5[m, :] = 0.0
            p6[m, :] = 0.0
            p7[m, :] = 0.0
            p8[m, :] = 0.0
            p9[m, :] = 0.0
            p10[m, :] = 0.0

            lti_min = []
            lti_max = []

            if NOISE_TON == 1:
                cp._default_memory_pool.free_all_blocks()
                tt_out = time.time()
                calc.function_farassat_1a(obj_par, obj_rot, temp_mics,
                                          "thickness")
                ttt_out = time.time()
                # print(f"computing time for thickness = "
                #       f"{ttt_out - tt_out}")

                tl_out = time.time()

                if TON_OPTIONS == 'f1a' or TON_OPTIONS == 'default':
                    calc.function_farassat_1a(obj_par, obj_rot, temp_mics,
                                              "loading")
                elif TON_OPTIONS == 'lowson':
                    calc.lowson(obj_par, obj_rot, temp_mics)

                tll_out = time.time()
                # print(f"computing time for loading = "
                #       f"{tll_out - tl_out}")

                calc.do_synthesis(obj_par, obj_rot, new_t2[m, :],
                                  ld_p2[m, :], tk_p2[m, :], tonal_p[m, :])

            if NOISE_BBN == 1:
                cp._default_memory_pool.free_all_blocks()
                geo.spherical_metrix(obj_par, obj_rot, temp_mics)
                calc.doppler_comp(obj_par, obj_rot, temp_mics)
                model.boundary_layer_thickness(obj_par, obj_rot)
                model.directivity(obj_rot)

                tb_out = time.time()

                if NOISE_TUR == 1:
                    p1[m, :], p2[m, :], p3[m, :] = model.bpm_tbl_te(
                        obj_par, obj_rot, obj_sig,
                        p1[m, :], p2[m, :], p3[m, :])
                if NOISE_LAM == 1:
                    p5[m, :] = model.bpm_lbl_vs(obj_par, obj_rot, obj_sig,
                                                p5[m, :])
                if NOISE_BLT == 1:
                    p6[m, :] = model.bpm_te_blt(obj_par, obj_rot, obj_sig,
                                                p6[m, :])
                if NOISE_TIP == 1:
                    p7[m, :] = model.bpm_tip_vs(obj_par, obj_rot, obj_sig,
                                                p7[m, :])
                if NOISE_BWI == 1:
                    p9[m, :], lti_min, lti_max = model.roger_bwi(obj_par, obj_rot, obj_sig,
                                                                 p9[m, :])

                p4[m, :] = p1[m, :] + p2[m, :] + p3[m, :]
                p8[m, :] = p4[m, :] + p5[m, :] + p6[m, :] + p7[m, :]
                p10[m, :] = p8[m, :] + p9[m, :]

                tbb_out = time.time()
                # print(f"computing time for BBN = "
                #       f"{tbb_out - tb_out}")

            for iii in range(len(lti_min)):
                test_f.write(f'{iii+1}\t{el[m]:.2f}\t{lti_min[iii]}\n'
                             f'{iii+1}\t{el[m]:.2f}\t{lti_max[iii]}\n')

            # print("dk2 =")
            # print(f'min = {dk2.min()}')
            # print(f'max = {dk2.max()}')
            # print('\n')
            # print("l2 =")
            # print(f'min = {l2.min()}')
            # print(f'max = {l2.max()}')
            # print('\n')
            # print("fcorr =")
            # print(f'min = {fcorr.min()}')
            # print(f'max = {fcorr.max()}')
            # print('\n')
            # print("fintg =")
            # print(f'min = {fintg.min()}')
            # print(f'max = {fintg.max()}')
            # print('\n')

        # te_noise = time.time()
        # print(f"computing time for noise calculation = "
        #       f"{te_noise - ts_noise}")

        obj_rot.new_t2 = new_t2
        obj_rot.ld_p2, obj_rot.tk_p2, obj_rot.tonal_p = ld_p2, tk_p2, tonal_p

        # Linearization for SPL-values.
        obj_rot.p1 = calc.calculation_spl(p1)
        obj_rot.p2 = calc.calculation_spl(p2)
        obj_rot.p3 = calc.calculation_spl(p3)
        obj_rot.p4 = calc.calculation_spl(p4)
        obj_rot.p5 = calc.calculation_spl(p5)
        obj_rot.p6 = calc.calculation_spl(p6)
        obj_rot.p7 = calc.calculation_spl(p7)
        obj_rot.p8 = calc.calculation_spl(p8)
        obj_rot.p9 = calc.calculation_spl(p9)
        obj_rot.p10 = calc.calculation_spl(p10)

        test_f.close()

        # Output files.
        # ts_out = time.time()
        out.write_tonal_noise(obj_par, obj_rot, obj_mic)
        out.write_broadband_noise(obj_par, obj_rot, obj_mic, obj_sig)
        # te_out = time.time()
        # print(f"computing time for Output = "
        #       f"{te_out - ts_out}")

    return None
