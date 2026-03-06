# import for system.
import os
import sys
# import for numerical calculation.
import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
# import for methods and functions.
import geometry as geo


def check_file(logs, inp_fn, arg_ignore):
    # Check whether the file exists.
    if not os.path.exists(inp_fn):
        logs.error('Error: There is no file "%s",'
                   ' please check the file.' % inp_fn)
        if arg_ignore:  # Ingnore the error.
            logs.error('But, now this error is ignored.')
        else:  # Stop the code.
            sys.exit('')

    return None


def get_configuration_parameters(args, logs, obj_par):
    # Initialize variables.
    par_right_name = {}
    par_wrong_name = {}
    par_wrong_value = {}

    # Check whether the file exists.
    check_file(logs, args.cfg, False)

    # Open the file.
    with open(args.cfg, 'r') as inp_f:
        # Read the file.
        for line in inp_f:
            # Check whether the line should be ignored.
            # when the line is empty or cmmented out.
            if line.startswith('#') or \
                    line.startswith(' ') or \
                    line == '\n':
                continue
            else:
                par_name = line.split('=')[0].strip()  # Parameter name
                par_value = line.split('=')[1].strip()  # Parameter value
                # Check whether the input parameter name is correct
                # and matched to one of the pre-defined paramteres.
                if hasattr(obj_par, par_name):
                    par_right_name[par_name] = par_value
                else:
                    par_wrong_name[par_name] = par_value
        # Close the file.

    # Check whether any variable is not correctly read.
    if (len(par_wrong_name) != 0):
        logs.error('Error: The following parameter names are not defined '
                   'correctly, please check the configuration file.')
        for name, value in par_wrong_name.items():
            logs.error('   - %s' % name)

        if args.ignore_cfg:
            logs.error('But, now this error is ignored.')
        else:
            sys.exit('')

    # Convert and set parameters
    type_hitnts = obj_par.__class__.__annotations__
    for par_name, par_value in par_right_name.items():
        try:
            expected_type = type_hitnts[par_name]
            if expected_type == int:
                converted_value = int(par_value)
            elif expected_type == float:
                converted_value = float(par_value)
            elif expected_type == str:
                converted_value = par_value

            setattr(obj_par, par_name, converted_value)

        except ValueError:
            par_wrong_value[par_name] = par_value
            continue

    if par_wrong_value:
        for name, value in par_wrong_value.items():
            logs.error('   - %s: %s' % (name, value))
        if args.ignore_cfg:
            logs.error('But, now this error is ignored.')
        else:
            sys.exit('')

    return None


def get_microphone_coordinates(obj_par, obj_mic):
    nmics = obj_par.nmics
    fname = os.path.join(obj_par.dir_inp, obj_par.MIC_FILENAME)
    datasets = np.genfromtxt(fname, skip_header=1,
                             delimiter=None, dtype=float)

    if nmics == 1:
        mic_x = datasets[0]
        mic_y = datasets[1]
        mic_z = datasets[2]
        mics = np.array([[mic_x, mic_y, mic_z]])
    else:
        mic_x = datasets[:, 0]
        mic_y = datasets[:, 1]
        mic_z = datasets[:, 2]
        mics = np.stack((mic_x, mic_y, mic_z), axis=1)

    # Convert to cupy ndarray for GPU computing.
    obj_mic.mics = cp.asarray(mics)

    return None


def get_deform(obj_par, obj_rot):
    NSTEP = obj_par.NSTEP
    NBLADE = obj_par.NBLADE
    NODEM, NODEN = obj_par.NODEM, obj_par.NODEN

    node = np.zeros((NSTEP, NBLADE, NODEN, NODEM, 3))

    fname = os.path.join(obj_par.dir_inp, obj_par.DEFORM_FILENAME)

    f = open(fname, 'r')
    f.readline()
    f.readline()

    # Import for blade geometry datas.
    for j in range(int(NODEN)):
        for i in range(int(NODEM)):
            def_temp = f.readline()
            node[:, :, j, i, 0] = float(def_temp.split()[0].strip())
            node[:, :, j, i, 1] = float(def_temp.split()[1].strip())
            node[:, :, j, i, 2] = float(def_temp.split()[2].strip())
    f.close()

    geo.apply_rotation_and_get_controlp(obj_par, obj_rot, node)

    return None


def get_loading_lowson(obj_par, obj_rot):
    nts = obj_par.nts
    NBLADE = obj_par.NBLADE

    NODEN = obj_par.NODEN - 1

    # Reading Point sources datasets.
    fname = os.path.join(obj_par.dir_inp, obj_par.LOADING_POINT_FILENAME)

    lowson_segs = np.zeros((nts, NBLADE, NODEN, 6))

    with open(fname, 'r') as f:
        next(f)
        for it in range(nts):
            next(f)
            for k in range(NBLADE):
                for j in range(NODEN):
                    lowson_temp = next(f).split()
                    lowson_segs[it, k, j, :] = np.array(
                        lowson_temp, dtype=float)

    obj_rot.lowson_segs = cp.asarray(lowson_segs)

    return None


def get_loading_f1a(obj_par, obj_rot):
    nts = obj_par.nts
    NBLADE = obj_par.NBLADE

    ld_nseg = obj_rot.ld_nseg

    fname = os.path.join(obj_par.dir_inp, obj_par.LOADING_FILENAME)

    cen_segs = np.zeros((nts, NBLADE, ld_nseg, 7))

    with open(fname, 'r') as f:
        next(f)
        for it in range(nts):
            next(f)
            for k in range(NBLADE):
                for j in range(ld_nseg):
                    loading_temp = next(f).split()
                    cen_segs[it, k, j, :] = np.array(loading_temp, dtype=float)

    obj_rot.cen_segs = cp.asarray(cen_segs)

    return None


# def get_tonal_datasets(obj_par, obj_rot):
#     nts = obj_par.nts
#     NBLADE = obj_par.NBLADE
#     NODEM, NODEN = obj_par.NODEM, obj_par.NODEN

#     ld_nseg = obj_rot.ld_nseg

#     # Reading Loading noise datasets ==========================================
#     fname1 = os.path.join(obj_par.dir_inp, obj_par.LOADING_FILENAME)

#     cen_segs = np.zeros((nts, NBLADE, ld_nseg, 7))

#     with open(fname1, 'r') as f1:
#         next(f1)
#         for it in range(nts):
#             next(f1)
#             for k in range(NBLADE):
#                 for j in range(ld_nseg):
#                     loading_temp = next(f1).split()
#                     cen_segs[it, k, j, :] = np.array(loading_temp, dtype=float)

#     obj_rot.cen_segs = cp.asarray(cen_segs)
#     # =========================================================================

#     # Reading Deformation datasets ============================================
#     NSTEP = obj_par.NSTEP

#     node = np.zeros((NSTEP, NBLADE, NODEN, NODEM, 3))

#     fname2 = os.path.join(obj_par.dir_inp, obj_par.DEFORM_FILENAME)

#     f2 = open(fname2, 'r')
#     f2.readline()
#     f2.readline()

#     # Import for blade geometry datas.
#     for j in range(int(NODEN)):
#         for i in range(int(NODEM)):
#             def_temp = f2.readline()
#             node[:, :, j, i, 0] = float(def_temp.split()[0].strip())
#             node[:, :, j, i, 1] = float(def_temp.split()[1].strip())
#             node[:, :, j, i, 2] = float(def_temp.split()[2].strip())
#     f2.close()

#     geo.apply_rotation_and_get_controlp(obj_par, obj_rot, node)
#     # =========================================================================

#     return None


def get_broadband(obj_par, obj_rot):
    nts = obj_par.nts
    NBLADE = obj_par.NBLADE
    NSEG = obj_par.NSEG

    x_le, y_le, z_le = obj_rot.x_le, obj_rot.y_le, obj_rot.z_le
    x_te, y_te, z_te = obj_rot.x_te, obj_rot.y_te, obj_rot.z_te
    C, L, alpha_geo = obj_rot.C, obj_rot.L, obj_rot.alpha_geo
    ufree = obj_rot.ufree
    alpha_star, alpha_zero = obj_rot.alpha_star, obj_rot.alpha_zero

    fname1 = os.path.join(obj_par.dir_inp, obj_par.BBN_CONDITION_FILENAME)

    f1 = open(fname1, 'r')
    next(f1)
    for it in range(nts+1):
        for k in range(NBLADE):
            next(f1)
            for j in range(NSEG):
                temp = f1.readline()
                x_le[it, k, j] = float(temp.split()[0].strip())
                y_le[it, k, j] = float(temp.split()[1].strip())
                z_le[it, k, j] = float(temp.split()[2].strip())
                x_te[it, k, j] = float(temp.split()[3].strip())
                y_te[it, k, j] = float(temp.split()[4].strip())
                z_te[it, k, j] = float(temp.split()[5].strip())
                C[it, k, j] = float(temp.split()[6].strip())
                L[it, k, j] = float(temp.split()[7].strip())
                alpha_geo[it, k, j] = float(temp.split()[8].strip())
                ufree[it, k, j] = float(temp.split()[9].strip())
                alpha_star[it, k, j] = float(temp.split()[10].strip())
                alpha_zero[it, k, j] = float(temp.split()[11].strip())
    f1.close()

    obj_rot.x_le = x_le
    obj_rot.y_le = y_le
    obj_rot.z_le = z_le
    obj_rot.x_te = x_te
    obj_rot.y_te = y_te
    obj_rot.z_te = z_te
    obj_rot.C = C
    obj_rot.L = L
    obj_rot.alpha_geo = alpha_geo
    obj_rot.ufree = ufree
    obj_rot.alpha_star = alpha_star
    obj_rot.alpha_zero = alpha_zero

    # Reading Turbulence characteristics parameters ===========================
    turb_leng, turb_ubar = obj_rot.turb_leng, obj_rot.turb_ubar
    turb_leng = cp.asnumpy(turb_leng)
    turb_ubar = cp.asnumpy(turb_ubar)

    fname2 = os.path.join(obj_par.dir_inp, obj_par.TURB_FILENAME)
    with open(fname2, 'r') as f2:
        for j in range(NSEG):
            temp = f2.readline()
            turb_leng[j] = float(temp.split()[0].strip())
            turb_ubar[j] = float(temp.split()[1].strip())

    turb_leng = cp.asarray(turb_leng)
    turb_ubar = cp.asarray(turb_ubar)

    obj_rot.turb_leng, obj_rot.turb_ubar = turb_leng, turb_ubar
    # =========================================================================

    return None
