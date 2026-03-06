# import for system.
import os
# import for numerical calculation.
import math
import cupy as cp
import numpy as np


class Parameter():
    # Program mode.
    MODE: str
    # Input/Output information.
    MIC_FILENAME: str
    DEFORM_FILENAME: str
    LOADING_FILENAME: str
    LOADING_POINT_FILENAME: str
    BBN_CONDITION_FILENAME: str
    TURB_FILENAME: str
    # FFT_FILENAME: str

    NOISE_TON: int
    TON_OPTIONS: str
    NOISE_BBN: int

    FRQ_BAND: str

    # Toggle options for B.B.N. calculation.
    BL_TRIP: int
    TE_RCORR: int
    NOISE_LAM: int
    NOISE_TUR: int
    NOISE_BLT: int
    NOISE_TIP: int
    NOISE_BWI: int

    # Blade condition definition.
    NSEG: int
    ASEG: int
    NI: int
    NJ: int
    NODEM: int
    NODEN: int
    NSTEP: int
    NPERI: int
    NBLADE: int
    R: float
    AOA_TIP: float
    TE_THICK: float
    TE_BLT_ANG: float

    # Flight condition definition.
    RPM: float
    MTIP: float
    AD: float

    # Microphone condition definition.
    OM: int
    ON: int
    ELV: float
    AZI: float

    # Flow condition definition
    MU: float
    RHO: float
    C0: float

    # Numerical methods definition.
    INTERP_LVL: int
    NINTERP: int

    def __init__(self):
        # Get a directory path for module operations.
        path_current = os.path.abspath(__file__)
        dir_current = os.path.dirname(path_current)

        self.dir_inp = os.path.join(dir_current, 'input')
        self.dir_out = os.path.join(dir_current, 'output')
        self.fft_inp = os.path.join(dir_current, 'FFT_input')
        self.fft_out = os.path.join(dir_current, 'FFT_output')

        self.MODE = 'default'
        self.MIC_FILENAME = ''
        self.DEFORM_FILENAME = ''
        self.LOADING_FILENAME = ''
        self.LOADING_POINT_FILENAME = ''
        self.BBN_CONDITION_FILENAME = ''
        self.TURB_FILENAME = ''
        # self.FFT_FILENAME = ''

        # Initial input parameters.
        self.NOISE_TON = 0
        self.TON_OPTIONS = ''
        self.NOISE_BBN = 0
        self.FRQ_BAND = ''
        self.BL_TRIP = 0
        self.TE_RCORR = 0
        self.NOISE_LAM = 0
        self.NOISE_TUR = 0
        self.NOISE_BLT = 0
        self.NOISE_TIP = 0
        self.NOISE_BWI = 0
        self.NSEG = 0
        self.ASEG = 0
        self.NI = 0
        self.NJ = 0
        self.NODEM = 0
        self.NODEN = 0
        self.NSTEP = 0
        self.NPERI = 0
        self.NBLADE = 0
        self.R = 0.0
        self.AOA_TIP = 0.0
        self.TE_THICK = 0.0
        self.TE_BLT_ANG = 0.0
        self.RPM = 0.0
        self.MTIP = 0.0
        self.AD = 0.0
        self.OM = 0
        self.ON = 0
        self.PHI = 0.0
        self.THETA = 0.0
        self.MU = 0.0
        self.RHO = 0.0
        self.C0 = 0.0
        self.INTERP_LVL = 4
        self.NINTERP = 0

        # Calculated values from input parameters.
        self.nmics = 0
        self.nfreq = 0
        self.rot_direc = 0

        self.nts = 0  # Number of time step.
        self.rps = 0.0  # Revolution per second.
        self.dt = 0.0

        self.mach_fwd = 0.0

        # self.omega = 0.0
        # self.dangle = 0.0
        # self.omega_dot = 0.0

    def set_variables(self, obj_sig):
        self.ASEG = (self.NI-1) * (self.NJ-1) * self.NBLADE
        self.MTIP = (self.RPM/60 * 2 * math.pi * self.R) / self.C0

        self.nmics = self.OM * self.ON
        if self.FRQ_BAND == 'one':
            self.nfreq = len(obj_sig.cen_ob)
        elif self.FRQ_BAND == 'third':
            self.nfreq = len(obj_sig.cen_ob3)
        elif self.FRQ_BAND == 'twelfth':
            self.nfreq = len(obj_sig.cen_ob12)

        if self.RPM >= 0:
            self.rot_direc = 1
        else:
            self.rot_direc = -1

        self.mach_fwd = self.MTIP * self.AD
        self.omega = float(self.MTIP * self.C0 / self.R)

        self.nts = int(self.NSTEP * self.NPERI)
        self.rps = float(self.RPM * self.rot_direc / 60)
        self.dt = 1 / (self.rps * self.NSTEP)

        # self.domega = float(360 / self.NSTEP)
        # self.omega_dot = self.RPM / 60 * 2 * cp.pi * self.rot_direc

        return None


class Signal():
    def __init__(self):
        # Frequency data for OB
        self.lbn_ob = [11, 22, 44, 88, 177, 355, 710, 1420, 2840, 5680, 11360]
        self.cen_ob = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000,
                       10000, 16000]
        self.ubn_ob = [22, 44, 88, 177, 355, 710, 1420, 2840, 5680, 11360,
                       22720]
        # Frequency data for 1/3 OB
        self.lbn_ob3 = [17.8, 22.4, 28.2, 35.5, 44.7, 56.2, 70.8, 89.1, 112.0,
                        141.0, 178.0, 224.0, 282.0, 355.0, 447.0, 562.0, 708.0,
                        891.0, 1120.0, 1410.0, 1780.0, 2240.0, 2820.0, 3550.0,
                        4470.0, 5620.0, 7080.0, 8910.0, 11200.0, 14130.0,
                        17780.0]
        self.cen_ob3 = [20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0,
                        160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0,
                        1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0,
                        5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0,
                        20000.0]
        self.ubn_ob3 = [22.4, 28.2, 35.5, 44.7, 56.2, 70.8, 89.1, 112.0, 141.0,
                        178.0, 224.0, 282.0, 355.0, 447.0, 562.0, 708.0, 891.0,
                        1120.0, 1410.0, 1780.0, 2240.0, 2820.0, 3550.0, 4470.0,
                        5620.0, 7080.0, 8910.0, 11200.0, 14130.0, 17780.0,
                        22390.0]
        # Frenquency data for 1/12 OB
        self.lbn_ob12 = [9.7, 10.3, 10.9, 11.5, 12.1, 12.8, 13.6, 14.6, 15.5,
                         16.5, 17.5, 18.5, 19.4, 20.6, 21.8, 22.9, 24.3, 25.7,
                         27.2, 29.1, 30.6, 32.5, 34.5, 36.4, 38.9, 41.3, 43.7,
                         46.2, 48.6, 51.5, 54.4, 58.3, 61.2, 65.1, 69.0, 72.9,
                         77.7, 82.6, 87.4, 92.3, 97.0, 103.0, 109.0, 115.0,
                         121.0, 128.0, 136.0, 146.0, 155.0, 165.0, 175.0,
                         185.0, 194.0, 206.0, 218.0, 229.0, 243.0, 257.0,
                         272.0, 291.0, 306.0, 325.0, 345.0, 364.0, 389.0,
                         413.0, 437.0, 462.0, 486.0, 515.0, 544.0, 583.0,
                         612.0, 651.0, 690.0, 729.0, 777.0, 826.0, 874.0,
                         923.0, 972.0, 1030.0, 1088.0, 1147.0, 1215.0, 1283.0,
                         1360.0, 1457.0, 1555.0, 1652.0, 1749.0, 1846.0,
                         1943.0, 2060.0, 2176.0, 2293.0, 2429.0, 2575.0,
                         2721.0, 2915.0, 3061.0, 3255.0, 3449.0, 3644.0,
                         3887.0, 4129.0, 4372.0, 4615.0, 4858.0, 5150.0,
                         5441.0, 5830.0, 6121.0, 6510.0, 6899.0, 7287.0,
                         7773.0, 8259.0, 8745.0, 9230.0, 9716.0, 10299.0,
                         10882.0, 11465.0, 12145.0, 12825.0, 13603.0, 14574.0,
                         15546.0, 16518.0, 17489.0, 18461.0, 19433.0, 21764.0,
                         24291.0, 27206.0]
        self.cen_ob12 = [10.0, 10.6, 11.2, 11.8, 12.5, 13.2, 14.0, 15.0, 16.0,
                         17.0, 18.0, 19.0, 20.0, 21.2, 22.4, 23.6, 25.0, 26.5,
                         28, 30.0, 31.5, 33.5, 35.5, 37.5, 40.0, 42.5, 45.0,
                         47.5, 50.0, 53.0, 56.0, 60.0, 63.0, 67.0, 71.0, 75.0,
                         80.0, 85.0, 90.0, 95.0, 100.0, 106.0, 112.0, 118.0,
                         125.0, 132.0, 140.0, 150.0, 160.0, 170.0, 180.0,
                         190.0, 200.0, 212.0, 224.0, 236.0, 250.0, 265.0,
                         280.0, 300.0, 315.0, 335.0, 355.0, 375.0, 400.0,
                         425.0, 450.0, 475.0, 500.0, 530.0, 560.0, 600.0,
                         630.0, 670.0, 710.0, 750.0, 800.0, 850.0, 900.0,
                         950.0, 1000.0, 1060.0, 1120.0, 1180.0, 1250.0, 1320.0,
                         1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0,
                         2000.0, 2120.0, 2240.0, 2360.0, 2500.0, 2650.0,
                         2800.0, 3000.0, 3150.0, 3350.0, 3550.0, 3750.0,
                         4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5300.0,
                         5600.0, 6000.0, 6300.0, 6700.0, 7100.0, 7500.0,
                         8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 10600.0,
                         11200.0, 11800.0, 12500.0, 13200.0, 14000.0, 15000.0,
                         16000.0, 17000.0, 18000.0, 19000.0, 20000.0, 22400.0,
                         25000.0, 28000.0]
        self.ubn_ob12 = [10.3, 10.9, 11.5, 12.1, 12.8, 13.6, 14.6, 15.5, 16.5,
                         17.5, 18.5, 19.4, 20.6, 21.8, 22.9, 24.3, 25.7, 27.2,
                         29.1, 30.6, 32.5, 34.5, 36.4, 38.9, 41.3, 43.7, 46.2,
                         48.6, 51.5, 54.4, 58.3, 61.2, 65.1, 69.0, 72.9, 77.7,
                         82.6, 87.4, 92.3, 97.0, 103.0, 109.0, 115.0, 121.0,
                         128.0, 136.0, 146.0, 155.0, 165.0, 175.0, 185.0,
                         194.0, 206.0, 218.0, 229.0, 243.0, 257.0, 272.0,
                         291.0, 306.0, 325.0, 345.0, 364.0, 389.0, 413.0,
                         437.0, 462.0, 486.0, 515.0, 544.0, 583.0, 612.0,
                         651.0, 690.0, 729.0, 777.0, 826.0, 874.0, 923.0,
                         972.0, 1030.0, 1088.0, 1147.0, 1215.0, 1283.0, 1360.0,
                         1457.0, 1555.0, 1652.0, 1749.0, 1846.0, 1943.0,
                         2060.0, 2176.0, 2293.0, 2429.0, 2575.0, 2721.0,
                         2915.0, 3061.0, 3255.0, 3449.0, 3644.0, 3887.0,
                         4129.0, 4372.0, 4615.0, 4858.0, 5150.0, 5441.0,
                         5830.0, 6121.0, 6510.0, 6899.0, 7287.0, 7773.0,
                         8259.0, 8745.0, 9230.0, 9716.0, 10299.0, 10882.0,
                         11465.0, 12145.0, 12825.0, 13603.0, 14574.0, 15546.0,
                         16518.0, 17489.0, 18461.0, 19433.0, 21764.0, 24291.0,
                         27206.0, 30606.0]
        # Set the frequency data as numpy arrays.
        self.cen_ob = cp.array(self.cen_ob)
        self.lbn_ob = cp.array(self.lbn_ob)
        self.ubn_ob = cp.array(self.ubn_ob)
        self.cen_ob3 = cp.array(self.cen_ob3)
        self.lbn_ob3 = cp.array(self.lbn_ob3)
        self.ubn_ob3 = cp.array(self.ubn_ob3)
        self.cen_ob12 = cp.array(self.cen_ob12)
        self.lbn_ob12 = cp.array(self.lbn_ob12)
        self.ubn_ob12 = cp.array(self.ubn_ob12)

        # Parameters for signal processing
        self.sw_dpl = False  # Switch for Doppler effect
        self.fft_ob_type = 'third'  # OB type: 'one', 'third', 'twelve'
        # Window function [~]: 'hamming', 'hanning'
        self.fft_wf_type = 'hamming'
        self.fft_blk_size = 2**10  # block size [\]
        self.fft_blk_len = 0.5  # Block length [s]
        self.fft_blk_ovl = 0.5  # Overlap length [-]

    def set_variables(self, obj_par):
        self.sw_dpl = obj_par.sig_sw_dpl
        self.spr_loss = obj_par.sig_spr_loss
        self.fft_ob_type = obj_par.sig_ob_type
        self.fft_wf_type = obj_par.sig_wf_type
        self.fft_blk_size = obj_par.sig_blk_size
        self.fft_blk_len = obj_par.sig_blk_len
        self.fft_blk_ovl = obj_par.sig_blk_ovl

        return None


class Microphone():
    def __init__(self):

        self.mics = []

        self.temp_mics = []

    def set_variables(self, obj_par):
        self.temp_mics = cp.zeros(3)

        return None


class Rotor():
    def __init__(self):
        # variables for Tonal noise.
        self.ld_nseg = 0.0
        self.tk_nseg = 0.0

        self.cen_segs = []  # (nts, NBLADE, NSEG, x or y or z)
        self.deform_segs = []
        self.lowson_segs = []  # for point sources, lowson methods.

        self.ld_src_t = []
        self.ld_obs_t = []
        self.ld_new_t = []
        self.ld_p = []
        self.ld_p_new = []
        self.ld_p_sum = []

        self.ls_p_new = []

        self.tk_src_t = []
        self.tk_obs_t = []
        self.tk_new_t = []
        self.tk_p = []
        self.tk_p_new = []
        self.tk_p_sum = []

        self.dt_ip = 0.0

        # Final value for Tonal noise.
        self.new_t2 = []
        self.ld_p2 = []
        self.tk_p2 = []
        self.tonal_p = []

        # variables for Broadband noise.
        self.x_le, self.y_le, self.z_le = [], [], []
        self.x_te, self.y_te, self.z_te = [], [], []
        self.C = []
        self.L = []
        self.alpha_geo = []
        self.ufree = []
        self.alpha_star = []
        self.alpha_zero = []
        self.r_le, self.elv_le, self.azi_le = [], [], []
        self.r_le, self.elv_te, self.azi_te = [], [], []
        self.cos_ksi_le, self.f_f0_le = [], []
        self.cos_ksi_te, self.f_f0_te = [], []

        self.h, self.psi = [], []

        self.turb_len = []
        self.turb_ubar = []

        self.re_chord = []
        self.mach = []

        self.delta_p, self.dstr_p, self.dstr_s = [], [], []
        self.dbar_h, self.dbar_l = [], []

        # SPL datas for broadband noise.
        self.spl_p = []
        self.spl_s = []
        self.spl_alpha = []
        self.p1, self.p2, self.p3, self.p4 = [], [], [], []
        self.p5, self.p6, self.p7, self.p8 = [], [], [], []
        self.p9, self.p10 = [], []

    def set_variables(self, obj_par):
        # variables for Tonal noise.
        self.ld_nseg = int(obj_par.ASEG / obj_par.NBLADE)
        self.ld_new_t = cp.zeros((obj_par.NINTERP), dtype=float)
        self.ld_p_new = cp.zeros((obj_par.NINTERP, obj_par.NBLADE,
                                  self.ld_nseg))
        self.ld_p_sum = cp.copy(self.ld_new_t)

        self.ls_p_new = cp.zeros((obj_par.NINTERP, obj_par.NBLADE,
                                  obj_par.NODEN - 1))

        self.tk_nseg = int((obj_par.NODEM - 1) * (obj_par.NODEN - 1))
        self.tk_new_t = cp.zeros((obj_par.NINTERP), dtype=float)
        self.tk_p_new = cp.zeros((obj_par.NINTERP, obj_par.NBLADE,
                                  self.tk_nseg))
        self.tk_p_sum = cp.copy(self.tk_new_t)

        self.new_t2 = cp.zeros((obj_par.nmics, obj_par.NINTERP))
        self.ld_p2 = cp.zeros((obj_par.nmics, obj_par.NINTERP))
        self.tk_p2 = cp.zeros((obj_par.nmics, obj_par.NINTERP))
        self.tonal_p = cp.zeros((obj_par.nmics, obj_par.NINTERP))

        # variables for Broadband noise.
        self.x_le = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.y_le = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.z_le = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.x_te = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.y_te = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.z_te = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                             dtype=cp.float64)
        self.C = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                          dtype=cp.float64)
        self.L = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                          dtype=cp.float64)
        self.alpha_geo = cp.zeros(
            (obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG), dtype=cp.float64)
        self.ufree = cp.zeros((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                              dtype=cp.float64)
        self.alpha_star = cp.zeros(
            (obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG), dtype=cp.float64)
        self.alpha_zero = cp.zeros(
            (obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG), dtype=cp.float64)

        self.h = cp.full((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                         obj_par.TE_THICK, dtype=cp.float64)
        self.psi = cp.full((obj_par.nts+1, obj_par.NBLADE, obj_par.NSEG),
                           obj_par.TE_BLT_ANG, dtype=cp.float64)

        self.turb_leng = cp.zeros(obj_par.NSEG)
        self.turb_ubar = cp.zeros(obj_par.NSEG)

        self.p1 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p2 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p3 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p4 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p5 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p6 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p7 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p8 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p9 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)
        self.p10 = cp.zeros((obj_par.nmics, obj_par.nfreq), dtype=cp.float64)

        return None
