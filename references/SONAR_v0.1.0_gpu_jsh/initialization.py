# import for system
import sys
from argparse import ArgumentParser
import logging


def explain_information(logs):
    logs.info(' ')
    logs.info('##############################################################')
    logs.info('#                                                            #')
    logs.info('#  NARC (Noise Analysis for Rotorcrafts Code)                #')
    logs.info('#  Version 0.1.0                                             #')
    logs.info('#                                                            #')
    logs.info('#  ADOLAB (Aeroacoustics & Design Optimization Lab.),        #')
    logs.info('#  Hanseo University, Korea.                                 #')
    logs.info('#                                                            #')
    logs.info('#  Main contributors                                         #')
    logs.info('#  Dr. Y. Jo (2023-current)                                  #')
    logs.info('#  Mr. S. Jeong (2023-current)                               #')
    logs.info('#                                                            #')
    logs.info('##############################################################')
    logs.info(' ')

    return None


def parse_arguments():
    """


    Returns
    -------
    class
        Arguments object.

    """
    # Take the input arguments using parser.
    parser = ArgumentParser()
    parser.add_argument("-f", "--file",
                        dest="cfg",
                        default="user_input.cfg",
                        help="read configuration parameters from FILE",
                        metavar="FILE")
    parser.add_argument("-q", "--quiet",
                        action="store_false",
                        dest="verbose",
                        default=True,
                        help="don't print messages")
    parser.add_argument("-l", "--log",
                        dest="log",
                        help="write log messages into FILE",
                        metavar="FILE")
    parser.add_argument("-e", "--err",
                        dest="err",
                        help="write error messages into FILE",
                        metavar="FILE")
    parser.add_argument("-i", "--ignore_cfg",
                        action="store_true",
                        dest="ignore_cfg",
                        default=False,
                        help="ignore all errors during reading cfg file")

    # Stop this code when the configuration file doesn't exist.
    if not parser.parse_args().cfg:
        print('Error: There is no configuration file specified. '
              'Please run this code with "-f" option '
              'and a configuration file name.')
        sys.exit('')

    return parser.parse_args()


def set_loggers(args):
    """

    Set loggers for logging messages during the process.

    Parameters
    ----------
    args : class
        Arguments object.

    Returns
    -------
    logs : class
        Loggers object.

    """
    # Initialize a logger.
    logs = logging.getLogger()
    logs.setLevel(level=logging.DEBUG)

    # Clear existing handlers
    if logs.hasHandlers():
        logs.handlers.clear()

    # Set the logger for printing info messages in the console.
    if args.verbose:
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        logs.addHandler(sh)
    # Set the logger for writing info messages in the log file.
    if args.log:
        fh = logging.FileHandler(args.log, mode='w')
        fh.setLevel(logging.INFO)
        logs.addHandler(fh)
    # Set the logger for writing error messages in the log file.
    if args.err:
        fh = logging.FileHandler(args.err, mode='w')
        fh.setLevel(logging.ERROR)
        logs.addHandler(fh)
    # Set the logger for writing error messages in the default file.
    # (saved in the execution path).
    else:
        fh = logging.FileHandler('error.log', mode='w')
        fh.setLevel(logging.ERROR)
        logs.addHandler(fh)

    return logs


def check_input_values(logs, obj_par):
    logs.info(' ')
    logs.info(f'              PROGRAM MODE: << {obj_par.MODE} >>')
    logs.info('  -- input parameters')
    logs.info(f'   - Speed of Sound                       = {obj_par.C0}')
    logs.info(
        f'   - Tip speed Mach                       = {obj_par.MTIP:.3f}')
    logs.info(f'   - Advance ratio                        = {obj_par.AD}')
    logs.info(f'   - Radius                               = {obj_par.R}')
    logs.info(f'   - Azimutal step (integer)              = {obj_par.NSTEP}')
    logs.info(f'   - Blade Num.                           = {obj_par.NBLADE}')
    logs.info(f'   - Blade segments Num.                  = {obj_par.NSEG}')
    logs.info(f'   - Surface chord node Num .             = {obj_par.NODEM}')
    logs.info(f'   - Surface span node Num.               = {obj_par.NODEN}')
    logs.info(f'   - Period                               = {obj_par.NPERI}')
    logs.info(f'   - Observer index i Num. (Map index i)  = {obj_par.OM}')
    logs.info(f'   - Observer index j Num. (Map index j)  = {obj_par.ON}')
    logs.info(
        f'   - Clockwise: 1, CounterClockwie: -1    = {obj_par.rot_direc}')
    logs.info(' ')
    logs.info('  -- Activated Toggle options,')
    logs.info('     Toggle option is ON: 1, OFF: 0')
    logs.info(
        f'   - Trip: {obj_par.BL_TRIP}         Round: {obj_par.TE_RCORR}')
    logs.info(f'   - Laminar BL vortex shedding noise: {obj_par.NOISE_LAM}')
    logs.info(
        f'   - Turbulence BL at trailing edge noise: {obj_par.NOISE_TUR}')
    logs.info(f'   - Trailing edge bluntness noise: {obj_par.NOISE_BLT}')
    logs.info(f'   - Trailing edge tip noise: {obj_par.NOISE_TIP}')
    logs.info(f'   - Blade-Wake interaction noise: {obj_par.NOISE_BWI}')
    logs.info(' ')

    return None
