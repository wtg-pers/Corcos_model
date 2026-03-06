# import for system.
import os
import sys
import time
from sys import platform
# import for methods and functions.
import variable as var
import file_input as inp
import initialization as init


def do_process(i1, args, logs):
    """
    Parameters
    ----------
    i1 : int
        First index of the process.
    args : class
        Arguments object.
    logs : class
        Loggers object.

    Returns
    -------
    objs : list
        List of objects.
    """

    # Initialize and set variables.
    i2 = 0

    # Create objects.
    i2 += 1
    logs.info(' --- %i.%i    Creating objects.' % (i1, i2))
    obj_par = var.Parameter()  # Parameter
    obj_rot = var.Rotor()
    obj_mic = var.Microphone()
    obj_sig = var.Signal()

    # Read parameter values from the configuration file.
    i2 += 1
    logs.info(' --- %i.%i    Reading user-defined parameters.' % (i1, i2))

    get_configuration_parameters(args, logs, obj_par)

    # Delete previous output files.
    i2 += 1
    logs.info(' --- %i.%i    Deleting previous results.' % (i1, i2))
    delete_previous_outputs(obj_par)

    # Set variables for the objects.
    i2 += 1
    logs.info(' --- %i.%i    Setting variables,' % (i1, i2))
    set_variables(i1, i2, logs, obj_par, obj_rot, obj_mic, obj_sig)
    # Set the objects into an object array.
    objs = []
    objs.append(obj_par)
    objs.append(obj_rot)
    objs.append(obj_mic)
    objs.append(obj_sig)

    # Read the aerodynamic datas for noise prediction.
    i2 += 1
    logs.info(' --- %i.%i    Reading Aerodynamic datas,' % (i1, i2))
    get_aerodynamic_datas(i1, i2, logs, objs)

    # Checking the input datas.
    init.check_input_values(logs, obj_par)

    return objs


def get_configuration_parameters(args, logs, obj_par):
    """
    Get configuration parameters from the input file.

    Parameters
    ----------
    None.

    Returns
    -------
    None.

    """
    inp.get_configuration_parameters(args, logs, obj_par)

    return None


def delete_previous_outputs(obj_par):
    """
    Delete output files previously generated. (v1.0)

    Parameters
    ----------
    obj_par : class
        Parameter object.

    Returns
    -------
    None.
    """
    # Set deletion command for the current system platform.
    if platform == 'linux' or platform == 'linux2' or platform == 'darwin':
        del_cmd = 'rm -rf *'
    elif platform == 'win32':
        del_cmd = 'del /q *'

    # Run the deletion operation corresponding to the output directory,
    if os.path.isdir(obj_par.dir_out):  # when the directory exists.
        os.system('cd %s && %s' % (obj_par.dir_out, del_cmd))
    else:  # when the directory is newly created.
        os.system('mkdir %s' % obj_par.dir_out)

    return None


def set_variables(i1, i2, logs, obj_par, obj_rot, obj_mic, obj_sig):
    """
    Parameters
    ----------
    i1 : int
        First index of the process.
    i2 : int
        Second index of the process.
    logs : class
        Loggers object.
    objs : list
        Object sets.

    Returns
    -------
    None.
    """

    # Initialize variables.
    i3 = 0  # Third index of the process

    # Set variables,
    # for initial datas.
    i3 += 1
    logs.info('  -- %i.%i.%i    for initialization.' % (i1, i2, i3))
    obj_par.set_variables(obj_sig)

    # for rotor objects.
    i3 += 1
    logs.info('  -- %i.%i.%i    for rotor.' % (i1, i2, i3))
    obj_rot.set_variables(obj_par)

    # for microphone object.
    i3 += 1
    logs.info('  -- %i.%i.%i    for microphones.' % (i1, i2, i3))
    obj_mic.set_variables(obj_par)

    return None


def get_aerodynamic_datas(i1, i2, logs, objs):
    """
    Parameters
    ----------
    i1 : int
        First index of the process.
    i2 : int
        Second index of the process.
    logs : class
        Loggers object.
    objs : list
        Object sets.

    Returns
    -------
    None.
    """
    obj_par = objs[0]
    obj_rot = objs[1]
    obj_mic = objs[2]
    # obj_sig = objs[3]
    # MODE = obj_par.MODE
    NOISE_TON = obj_par.NOISE_TON
    TON_OPTIONS = obj_par.TON_OPTIONS
    NOISE_BBN = obj_par.NOISE_BBN

    # print(obj_par.ASEG)

    # Initialize variables.
    i3 = 0  # Third index of the process

    # Set variables,
    # for initial datas.
    i3 += 1
    logs.info('  -- %i.%i.%i    for microphone coordinates.' %
              (i1, i2, i3))
    ts_mic = time.time()
    inp.get_microphone_coordinates(obj_par, obj_mic)
    te_mic = time.time()

    # print(f"computing time for Input Microphone datas = "
    #       f"{te_mic - ts_mic}")

    if NOISE_TON == 1:
        i3 += 1
        logs.info('  -- %i.%i.%i    for Deformation.' % (i1, i2, i3))
        ts_deform = time.time()
        inp.get_deform(obj_par, obj_rot)
        te_deform = time.time()

        i3 += 1
        logs.info('  -- %i.%i.%i    for loading.' % (i1, i2, i3))
        ts_ld = time.time()
        if TON_OPTIONS == 'f1a' or TON_OPTIONS == 'default':
            inp.get_loading_f1a(obj_par, obj_rot)
        elif TON_OPTIONS == 'lowson':
            inp.get_loading_lowson(obj_par, obj_rot)
        te_ld = time.time()

        # print(f"computing time for Input Deform = "
        #       f"{te_deform - ts_deform}")
        # print(f"computing time for Input Loading = "
        #       f"{te_ld - ts_ld}")

    if NOISE_BBN == 1:
        i3 += 1
        logs.info('  -- %i.%i.%i    for broadband.' % (i1, i2, i3))
        ts_bbn = time.time()
        inp.get_broadband(obj_par, obj_rot)
        te_bbn = time.time()

        # print(f"computing time for Input BBN = "
        #       f"{te_bbn - ts_bbn}")

    return None
