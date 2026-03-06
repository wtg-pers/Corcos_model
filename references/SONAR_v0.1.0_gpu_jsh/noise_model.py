# import for system.
import os
import time
# import for numberical calculation.
import numpy as np
import cupy as cp
import scipy.special as sp
import math as m
from multiprocessing import Pool
from scipy.integrate import quad


def boundary_layer_thickness(obj_par, obj_rot):
    """
    Variables
    delta : Boundary layer thickness
    delta0 : Zero angle of attack delta
    dstr : Displacement thickness, delta^*
    dstr0 : Zero angle of attack dstr
    dstrs : Suction side dstr
    dstrp : Pressure side dstr
    """

    re_chord = obj_rot.re_chord
    C = obj_rot.C
    a_star, a_zero = obj_rot.alpha_star, obj_rot.alpha_zero
    alpha_star = a_star - a_zero

    BL_TRIP = obj_par.BL_TRIP

    # Compute zero angle of attack boundary layer thickness.
    # ------------------------------------------------------
    if ((BL_TRIP == 1) or (BL_TRIP == 2)):
        delta0 = 10.0**(1.8920 - 0.9045 * cp.log10(re_chord) +
                        0.0596 * cp.log10(re_chord)**2.0) * C

        if (BL_TRIP == 2):
            delta0 *= 0.6
    else:
        delta0 = 10.0**(1.6569 - 0.9045 * cp.log10(re_chord) +
                        0.0596 * cp.log10(re_chord)**2.0) * C

    # Compute pressure side BL thickness.
    # -----------------------------------
    delta_p = 10.0**(-0.04175 * alpha_star + 0.00106 * alpha_star**2) * delta0

    # Compute zero angle of aatack displacement thickness.
    # !!! Check the re_c criteria is 0.3x10^6, not 3.0x10^6 !!!
    # ---------------------------------------------------------
    if ((BL_TRIP == 1) or (BL_TRIP == 2)):
        dstr0 = cp.where(re_chord <= 0.3E6,
                         0.0601 * re_chord**(-0.114) * C,
                         10.0**(3.411 - 1.5397 *
                                cp.log10(re_chord) + 0.1059 *
                                cp.log10(re_chord)**2.0) * C
                         )

        if (BL_TRIP == 2):
            dstr0 *= 0.6

    else:
        dstr0 = 10.0**(3.0187 - 1.5397 * cp.log10(re_chord) +
                       0.1059 * cp.log10(re_chord)**2.0) * C

    # Pressure side displacement thickness.
    # -------------------------------------
    dstr_p = 10.0**(-0.0432 * alpha_star + 0.00113 * alpha_star**2.0) * dstr0

    if (BL_TRIP == 3):
        dstr_p *= 1.48

    # Suction side displacement thickness.
    # ------------------------------------
    if (BL_TRIP == 1):
        cond1 = alpha_star <= 5.0
        cond2 = (alpha_star > 5.0) & (alpha_star <= 12.5)
        # cond3 = alpha_star > 12.5
        dstr_s = cp.where(
            cond1,
            10.0**(0.0679 * alpha_star) * dstr0,
            cp.where(
                cond2,
                0.381 * 10.0**(0.1516 * alpha_star) * dstr0,
                14.296 * 10.0**(0.0258 * alpha_star) * dstr0
            )
        )

    else:
        cond1 = alpha_star <= 7.5
        cond2 = (alpha_star > 7.5) & (alpha_star <= 12.5)
        # cond3 = alpha_star > 12.5
        dstr_s = cp.where(
            cond1,
            10.0**(0.0679 * alpha_star) * dstr0,
            cp.where(
                cond2,
                0.0162 * 10.0**(0.3066 * alpha_star) * dstr0,
                52.420 * 10.0**(0.0258 * alpha_star) * dstr0
            )
        )

    obj_rot.delta_p = delta_p
    obj_rot.dstr_p = dstr_p
    obj_rot.dstr_s = dstr_s

    # print(delta_p)

    return None


def directivity(obj_rot):
    mach = obj_rot.mach
    azi_le, elv_le = obj_rot.azi_le, obj_rot.elv_le
    cos_ksi_le = obj_rot.cos_ksi_le

    azi_te, elv_te = obj_rot.azi_te, obj_rot.elv_te
    cos_ksi_te = obj_rot.cos_ksi_te

    azir_le = cp.deg2rad(azi_le)
    elvr_le = cp.deg2rad(elv_le)

    azir_te = cp.deg2rad(azi_te)
    elvr_te = cp.deg2rad(elv_te)

    # dbar_h = 2.0 * (cp.sin(azir / 2.0)**2 *
    #                 cp.sin(elvr)**2) / (1.0 - mach * cos_ksi)**4
    dbar_h_le = 2.0 * (cp.sin(elvr_le / 2.0)**2 *
                       cp.sin(azir_le)**2) / (1.0 - mach * cos_ksi_le)**4
    dbar_h_te = 2.0 * (cp.sin(elvr_te / 2.0)**2 *
                       cp.sin(azir_te)**2) / (1.0 - mach * cos_ksi_te)**4

    dbar_l_le = (cp.sin(azir_le)**2 * cp.sin(elvr_le)**2
                 ) / (1.0 - mach * cos_ksi_le)**4
    dbar_l_te = (cp.sin(azir_te)**2 * cp.sin(elvr_te)**2
                 ) / (1.0 - mach * cos_ksi_te)**4

    obj_rot.dbar_h_le = dbar_h_le
    obj_rot.dbar_l_le = dbar_l_le

    obj_rot.dbar_h_te = dbar_h_te
    obj_rot.dbar_l_te = dbar_l_te

    return None


def a0comp(re_chord):
    # This subroutine determines where the A-curve take on A value of -20 dB.
    # -----------------------------------------------------------------------
    # re_chord = obj_rot.re_chord

    cond1 = re_chord < 9.52E4
    cond2 = (re_chord >= 9.52E4) & (re_chord <= 8.57E5)
    cond3 = re_chord > 8.57E5

    a0 = cp.where(
        cond1,
        0.57,
        cp.where(
            cond2,
            (-9.57E-13) * (re_chord - 8.57E5)**2.0 + 1.13,
            1.13)
    )

    return a0


def amin(a):
    # This subroutine defines the curve fit corresponding to the A-curve
    # for the minimum allowed reynolds number.
    # ------------------------------------------------------------------
    x = abs(a)

    cond1 = x < 0.204
    cond2 = (x >= 0.204) & (x <= 0.244)
    cond3 = x > 0.244

    amina = cp.where(
        cond1,
        cp.sqrt(67.552 - 886.788 * x**2.0) - 8.219,
        cp.where(
            cond2,
            -32.665 * x + 3.981,
            -142.795 * x**3.0 + 103.656 * x**2.0 - 57.757 * x + 6.006
        )
    )

    return amina


def amax(a):
    # This subroutine defines the curve fit corresponding to the A-curve
    # for the maximum allowed reynolds number.
    # ------------------------------------------------------------------
    x = abs(a)

    cond1 = x < 0.13
    cond2 = (x >= 0.13) & (x <= 0.321)
    cond3 = x > 0.321

    amaxa = cp.where(
        cond1,
        cp.sqrt(67.552 - 886.788 * x**2.0) - 8.219,
        cp.where(
            cond2,
            -15.901 * x + 1.098,
            -4.669 * x**3.0 + 3.491 * x**2.0 - 16.699 * x + 1.149
        )
    )

    return amaxa


def bmin(b):
    # This subroutine defines the curve fit corresponding to the B-curve
    # for the minimum allowed reynolds number.
    # ------------------------------------------------------------------
    x = abs(b)

    cond1 = x < 0.13
    cond2 = (x >= 0.13) & (x <= 0.145)
    cond3 = x > 0.145

    bminb = cp.where(
        cond1,
        cp.sqrt(16.888 - 886.788 * x**2.0) - 4.109,
        cp.where(
            cond2,
            -83.607 * x + 8.138,
            -817.81 * x**3.0 + 355.21 * x**2.0 - 135.024 * x + 10.619
        )
    )

    return bminb


def bmax(b):
    # This subroutine defines the curve fit corresponding to the B-curve
    # for the maximum allowed reynolds number.
    # ------------------------------------------------------------------
    x = abs(b)

    cond1 = x < 0.1
    cond2 = (x >= 0.1) & (x <= 0.187)
    cond3 = x > 0.187

    bmaxb = cp.where(
        cond1,
        cp.sqrt(16.888 - 886.788 * x**2.0) - 4.109,
        cp.where(
            cond2,
            -31.33 * x + 1.854,
            -80.541 * x**3.0 + 44.174 * x**2.0 - 39.381 * x + 2.344
        )
    )

    return bmaxb


def g5comp(hdstar, eta):
    m = 0.0
    mu = 0.0

    cond1 = (hdstar < 0.25)
    cond2 = (hdstar >= 0.25) & (hdstar < 0.62)
    cond3 = (hdstar >= 0.62) & (hdstar < 1.15)
    mu = cp.where(cond1, 0.1221,
                  cp.where(cond2, -0.2175 * hdstar + 0.1755,
                           cp.where(cond3, -0.0308 * hdstar + 0.0596,
                                    0.0242)))

    cond1 = (hdstar <= 0.02)
    cond2 = (hdstar > 0.02) & (hdstar <= 0.5)
    cond3 = (hdstar > 0.5) & (hdstar <= 0.62)
    cond4 = (hdstar > 0.62) & (hdstar <= 1.15)
    cond5 = (hdstar > 1.15) & (hdstar < 1.2)
    m = cp.where(cond1, 0.0,
                 cp.where(cond2, 68.724 * hdstar - 1.35,
                          cp.where(cond3, 308.475 * hdstar - 121.23,
                                   cp.where(cond4, 224.811 * hdstar - 69.354,
                                            cp.where(cond5,
                                                     1583.28 *
                                                     hdstar - 1631.592,
                                                     268.344)))))

    cond = (m < 0.0)
    m = cp.where(cond, 0.0, m)

    eta0 = -cp.sqrt((m**2 * mu**4) / (6.25 + m**2 * mu**2))
    k = 2.5 * cp.sqrt(1.0 - (eta0 / mu)**2) - 2.5 - m * eta0

    cond1 = (eta < eta0)
    cond2 = (eta >= eta0) & (eta < 0.0)
    cond3 = (eta >= 0.0) & (eta < 0.03616)
    g5 = cp.where(cond1, m * eta + k,
                  cp.where(cond2, 2.5 * cp.sqrt(1.0 - (eta / mu)**2) - 2.5,
                           cp.where(cond3, cp.sqrt(
                               1.5625 - 1194.99 * eta**2) - 1.25,
                               -155.543 * eta + 4.375)))

    return g5


def bpm_tbl_te(obj_par, obj_rot, obj_sig, p1, p2, p3):
    NBLADE, NSEG = obj_par.NBLADE, obj_par.NSEG
    NPERI = obj_par.NPERI
    NSTEP = obj_par.NSTEP
    nts = obj_par.nts
    nmics = obj_par.nmics
    MU, RHO = obj_par.MU, obj_par.RHO
    RPM = obj_par.RPM

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    nfreq = obj_par.nfreq

    r = obj_rot.r_te
    L = obj_rot.L
    mach = obj_rot.mach
    re_chord = obj_rot.re_chord
    ufree = obj_rot.ufree
    a_star, a_zero = obj_rot.alpha_star, obj_rot.alpha_zero
    f_f0 = obj_rot.f_f0_te

    alpha_star = a_star - a_zero

    pgcomp = cp.sqrt(1.0 - mach**2)

    # Compute boundary layer thickness.
    dstr_p, dstr_s = obj_rot.dstr_p, obj_rot.dstr_s
    # ---------------------------------

    # Compute directivity function.
    dbar_h, dbar_l = obj_rot.dbar_h_te, obj_rot.dbar_l_te
    # -----------------------------

    # Calculate the Reynolds number based on pressure and
    # suction displacement thickness.
    rdstr_p = dstr_p * ufree / (MU / RHO)
    rdstr_s = dstr_s * ufree / (MU / RHO)
    # ---------------------------------------------------

    # Determine peak Strouhal number to be used for
    # 'A' and 'B' curve calculations.
    # ---------------------------------------------
    st1 = 0.02 * mach**(-0.6)

    cond1 = alpha_star < 1.333
    cond2 = (alpha_star >= 1.333) & (alpha_star <= 12.5)
    # cond3 = alpha_star > 12.5
    st2 = cp.where(cond1,
                   st1,
                   cp.where(cond2,
                            st1 * 10.0**(0.0054 * (alpha_star - 1.333)**2),
                            4.72 * st1
                            )
                   )

    st1prim = (st1 + st2) / 2.0

    a0 = a0comp(re_chord)
    a02 = a0comp(3.0 * re_chord)

    # Evaluate minimum and maximum 'A' curves at a0.
    # Shape is (nts, NBLADE, NSEG)
    # ----------------------------------------------
    amina0 = amin(a0)
    amaxa0 = amax(a0)
    amina02 = amin(a02)
    amaxa02 = amax(a02)

    # Compute 'A' max/min ratio.
    # --------------------------
    ara0 = (20.0 + amina0) / (amina0 - amaxa0)
    ara02 = (20.0 + amina02) / (amina02 - amaxa02)

    # Compute b0 to be used in 'B' curve calculation.
    # -----------------------------------------------
    cond1 = re_chord < 9.52e+4
    cond2 = (re_chord >= 9.52e+4) & (re_chord <= 8.57e+5)
    # cond3 = re_chord >= 8.57E5
    b0 = cp.where(cond1,
                  0.3,
                  cp.where(cond2,
                           (-4.48e-13) * (re_chord - 8.57e+5)**2.0 + 0.56,
                           0.56
                           )
                  )

    # Evaluate minimum and maximum 'B' curves at b0.
    # ----------------------------------------------
    bminb0 = bmin(b0)
    bmaxb0 = bmax(b0)

    # Compute 'B' max/min ratio.
    # --------------------------
    brb0 = (20.0 + bminb0) / (bminb0 - bmaxb0)

    # For each center frequency, compute an 'A' prediction
    # for the pressure side.
    # ----------------------------------------------------
    """
    id : jsh
    have been changed the value assigned to 'stpeak'
    st2 -> st1 ref: fortran SOURCE code(no rotating)
    """
    if (RPM == 0.0):
        stpeak = st1
    else:
        stpeak = st2
    # stpeak = st1  # FOR FIXED WING
    # stpeak = st2  # FOR ROTATIONG WING
    # stpeak = st1prim
    # stpeak shape is (nts, NBLADE, NSEG)
    du_temp = (dstr_p / ufree)
    st_p = cen_ob[:, cp.newaxis, cp.newaxis, cp.newaxis
                  ] * du_temp[cp.newaxis, :, :, :]
    a = cp.log10(st_p / stpeak[cp.newaxis, :, :, :])
    # st_p, a shape is (nfreq, nts, NBLADE, NSEG)

    amina = amin(a)
    amaxa = amax(a)
    aa = amina + ara0[cp.newaxis, :, :, :] * (amaxa - amina)

    cond1 = re_chord < 2.47e+5
    cond2 = (re_chord >= 2.47e+5) & (re_chord <= 8.0e+5)
    # cond3 = re_chord >= 8.0E5
    k1 = cp.where(cond1,
                  -4.31 * cp.log10(re_chord) + 156.3,
                  cp.where(cond2,
                           -9.0 * cp.log10(re_chord) + 181.6,
                           128.5
                           )
                  )

    cond1 = rdstr_p <= 5000.0
    # cond2 = rdstr_p > 5000.0
    delk1 = cp.where(cond1,
                     alpha_star * (1.43 * cp.log10(rdstr_p) - 5.29),
                     0.0)

    temp_p = k1 - 3 + 10 * cp.log10(dstr_p * mach**5 * dbar_h * L /
                                    r**2) + delk1
    spl_p = aa + temp_p[cp.newaxis, :, :, :]

    gamma = 27.094 * mach + 3.310
    gamma0 = 23.430 * mach + 4.651
    beta = 72.650 * mach + 10.740
    beta0 = -34.190 * mach - 13.820

    cond1 = alpha_star < (gamma0 - gamma)
    cond2 = (alpha_star >= (gamma0 - gamma)) & (alpha_star <= (gamma0 + gamma))
    # cond3 = alpha_star > (gamma0 + gamma)
    k2 = cp.where(
        cond1,
        -1000.0,
        cp.where(
            cond2,
            cp.sqrt(beta**2.0 - (beta / gamma)**2.0 *
                    (alpha_star - gamma0)**2.0) + beta0,
            12.0
        )
    )
    k2 += k1  # (nts, NBLADE, NSEG)

    du_temp = (dstr_s / ufree)
    st_s = cen_ob[:, cp.newaxis, cp.newaxis, cp.newaxis
                  ] * du_temp[cp.newaxis, :, :, :]

    # Check for 'A' computation for suction side.
    # -------------------------------------------
    xcheck = gamma0
    switch = False

    cond = ((alpha_star >= xcheck) | (alpha_star > 12.5))
    switch = cp.where(cond, True, False)

    if RPM == 0.0:
        a = cp.where(switch, a, cp.log10(st_s / st1prim[cp.newaxis, :, :, :]))
    else:
        a = cp.where(switch, a, cp.log10(st_s / stpeak[cp.newaxis, :, :, :]))

    amina = cp.where(switch, amina, amin(a))

    amaxa = cp.where(switch, amaxa, amax(a))

    AA = cp.where(switch,
                  aa,
                  amina + ara0[cp.newaxis, :, :, :] * (amaxa - amina)
                  )

    temp_s = k1 - 3 + 10 * cp.log10(dstr_s * mach**5 * dbar_h * L / r**2)

    spl_s = cp.where(switch,
                     10 * cp.log10(dstr_s * mach**5 * dbar_l * L / r**2),
                     AA + temp_s[cp.newaxis, :, :, :]
                     )

    spl_p = cp.where(switch,
                     10 * cp.log10(dstr_s * mach**5 * dbar_l * L / r**2),
                     spl_p)

    # 'B' curve computation.
    # ----------------------
    b = cp.where(switch,
                 abs(cp.log10(st_s / st2[cp.newaxis, :, :, :])),
                 cp.log10(st_s / st2[cp.newaxis, :, :, :]))

    minb = cp.where(switch,
                    amin(b),
                    bmin(b))
    maxb = cp.where(switch,
                    amax(b),
                    bmax(b))
    bb = cp.where(switch,
                  minb + ara02 * (maxb - minb),
                  minb + brb0 * (maxb - minb))

    spl_alpha = cp.where(switch,
                         bb + k2 + 10 * cp.log10(dstr_s * mach**5 * dbar_l *
                                                 L / r**2),
                         bb + k2 + 10 * cp.log10(dstr_s * mach**5 * dbar_h *
                                                 L / r**2))

    # Prantdl-Glauert corection
    # !!!CHECK For ROTATION WING ONLY!!!
    # ----------------------------------
    if (RPM != 0.0):
        spl_p /= pgcomp[cp.newaxis, :, :, :]**2
        spl_s /= pgcomp[cp.newaxis, :, :, :]**2
        spl_alpha /= pgcomp[cp.newaxis, :, :, :]**2

    # Sum all contributions from 'A' and 'B' on both
    # pressure and suction side on a mean-square pressure basis.
    # ----------------------------------------------------------
    condp = spl_p < -100.0
    spl_p = cp.where(
        condp,
        -100.0,
        spl_p)

    conds = spl_s < -100.0
    spl_s = cp.where(
        conds,
        -100.0,
        spl_s)

    condalpha = spl_alpha < -100.0
    spl_alpha = cp.where(
        condalpha,
        -100.0,
        spl_alpha)

    # Calculation results of p'
    # scaled_f_f0 = f_f0 / (nts+1)

    scaled_f_f0 = f_f0 / NSTEP
    scaled_f_f0 = scaled_f_f0[cp.newaxis, :, :, :]

    p1 = cp.sum(10.0**(spl_p / 10.0) * scaled_f_f0, axis=(3, 2, 1))
    p1 /= NPERI
    # p1[:] += result_p

    p2 = cp.sum(10.0**(spl_s / 10.0) * scaled_f_f0, axis=(3, 2, 1))
    p2 /= NPERI
    # p2[:] += result_s

    p3 = cp.sum(10.0**(spl_alpha / 10.0) * scaled_f_f0,
                axis=(3, 2, 1))
    p3 /= NPERI
    # p3[:] += result_alpha

    return p1, p2, p3


def bpm_lbl_vs(obj_par, obj_rot, obj_sig, p5):
    NSTEP = obj_par.NSTEP
    NPERI = obj_par.NPERI
    nts = obj_par.nts
    RPM = obj_par.RPM

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    nfreq = obj_par.nfreq

    r = obj_rot.r_te
    L = obj_rot.L
    mach = obj_rot.mach
    re_chord = obj_rot.re_chord
    ufree = obj_rot.ufree
    a_star, a_zero = obj_rot.alpha_star, obj_rot.alpha_zero
    f_f0 = obj_rot.f_f0_te

    # Compute the effective angle of attack.
    alpha_star = a_star - a_zero

    # Compute boundary layer thickness.
    delta_p = obj_rot.delta_p
    # ---------------------------------

    # Compute directivity function.
    dbar_h = obj_rot.dbar_h_te
    # -----------------------------

    # Compute reference strouhal number.
    # ----------------------------------
    cond1 = (re_chord <= 1.3E5)
    cond2 = (re_chord > 1.3E5) & (re_chord <= 4E5)
    # cond3 = (re_chord > 4E5)
    st1prim = cp.where(cond1, 0.18,
                       cp.where(cond2, 0.001756 * re_chord**0.3931, 0.28))

    stpkprm = st1prim * 10.0**(-0.04 * alpha_star)

    # Compute reference reynolds number.
    # ----------------------------------
    cond = (alpha_star <= 3.0)
    re_chord0 = cp.where(cond, 10.0**(0.215 * alpha_star + 4.978),
                         10.0**(0.120 * alpha_star + 5.263))

    # Compute peak scaled spectrum level.
    # -----------------------------------
    d = re_chord / re_chord0

    cond1 = (d <= 0.3237)
    cond2 = (d > 0.3237) & (d <= 0.5689)
    cond3 = (d > 0.5689) & (d <= 1.7579)
    cond4 = (d > 1.7579) & (d <= 3.0889)
    # cond5 = (d > 3.0889)
    g2 = cp.where(cond1, 77.852 * cp.log10(d) + 15.328,
                  cp.where(cond2, 65.188 * cp.log10(d) + 9.125,
                           cp.where(cond3, -144.052 * cp.log10(d)**2.0,
                                    cp.where(cond4,
                                             -65.188 * cp.log10(d) + 9.125,
                                             -77.852 * cp.log10(d) + 15.328))))

    g3 = 171.04 - 3.03 * alpha_star

    scale = 10.0 * cp.log10(delta_p * mach**5 * dbar_h * L / r**2.0)

    # Compute scaled SPL for each strouhal number.
    # --------------------------------------------
    du_temp = (delta_p / ufree)
    stprim = cen_ob[:, cp.newaxis, cp.newaxis, cp.newaxis
                    ] * du_temp[cp.newaxis, :, :, :]
    e = stprim / stpkprm[cp.newaxis, :, :, :]

    cond1 = (e < 0.5974)
    cond2 = (e >= 0.5974) & (e < 0.8545)
    cond3 = (e >= 0.8545) & (e < 1.17)
    cond4 = (e >= 1.17) & (e < 1.674)
    # cond5 = (e > 1.674)
    g1 = cp.where(cond1, 39.80 * cp.log10(e) - 11.12,
                  cp.where(cond2, 98.409 * cp.log10(e) + 2.0,
                           cp.where(cond3, -5.076 +
                                    cp.sqrt(2.484 - 506.25 * cp.log10(e)**2.0),
                                    cp.where(cond4,
                                             -98.409 * cp.log10(e) + 2.0,
                                             -39.80 * cp.log10(e) - 11.12))))

    spl_lbl = g1 + g2[cp.newaxis, :, :, :] + g3[cp.newaxis, :, :, :] + \
        scale[cp.newaxis, :, :, :]

    # Calculation results of p'
    scaled_f_f0 = f_f0 / NSTEP
    scaled_f_f0 = scaled_f_f0[cp.newaxis, :, :, :]

    p5 = cp.sum(10**(spl_lbl / 10) * scaled_f_f0, axis=(1, 2, 3))
    p5 /= NPERI

    return p5


def bpm_te_blt(obj_par, obj_rot, obj_sig, p6):
    NBLADE, NSEG = obj_par.NBLADE, obj_par.NSEG
    NSTEP = obj_par.NSTEP
    NPERI = obj_par.NPERI
    nts = obj_par.nts
    RPM = obj_par.RPM

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    nfreq = obj_par.nfreq

    r = obj_rot.r_te
    L = obj_rot.L
    h, psi = obj_rot.h, obj_rot.psi
    mach = obj_rot.mach
    ufree = obj_rot.ufree

    f_f0 = obj_rot.f_f0_te

    f4temp = 0.0
    g4 = 0.0
    g5 = 0.0
    g50 = 0.0
    g514 = 0.0

    # Compute boundary layer thickness.
    dstr_p, dstr_s = obj_rot.dstr_p, obj_rot.dstr_s
    # ---------------------------------

    # Compute directivity function.
    dbar_h = obj_rot.dbar_h_te
    # -----------------------------

    # Compute averaged displacement thickness.
    dstr_avg = (dstr_p + dstr_s) / 2
    hdstar = h / dstr_avg
    dstarh = 1.0 / hdstar
    # ----------------------------------------

    # Compute peak Strouhal number.
    # -----------------------------
    aterm = 0.212 - 0.0045 * psi

    cond = (hdstar >= 0.2)
    stpeak = cp.where(cond,
                      aterm / (1.0 + 0.235 * dstarh - 0.0132 * dstarh**2),
                      0.1 * hdstar + 0.095 - 0.00243 * psi)

    # Compute scaled spectrum level.
    # ------------------------------
    cond = (hdstar <= 5.0)
    g4 = cp.where(cond, 17.5 * cp.log10(hdstar) + 157.5 - 1.114 * psi,
                  169.7 - 1.114 * psi)

    # For each frequency, compute spectrum shape referenced to 0 dB.
    # --------------------------------------------------------------
    hu_temp = h / ufree
    stppp = cen_ob[:, cp.newaxis, cp.newaxis, cp.newaxis
                   ] * hu_temp[cp.newaxis, :, :, :]
    eta = cp.log10(stppp / stpeak[cp.newaxis, :, :, :])

    hdstarl = hdstar
    # hdstar_l = hdstar_l / 7.0

    g514 = g5comp(hdstarl[cp.newaxis, :, :, :], eta)

    hdstar_p = 6.724 * hdstar**2 - 4.019 * hdstar + 1.107
    # hdstar_p = hdstar_p / 7.0

    g50 = g5comp(hdstar_p[cp.newaxis, :, :, :], eta)

    g5 = g50 + 0.0714 * psi[cp.newaxis, :, :, :] * (g514 - g50)

    cond = (g5 > 0.0)
    g5 = cp.where(cond, 0.0, g5)

    hdstar025 = cp.full((nfreq, (nts+1), NBLADE, NSEG), 0.25)
    f4temp = g5comp(hdstar025, eta)

    cond = (g5 > f4temp)
    g5 = cp.where(cond, f4temp, g5)

    scale = 10.0 * cp.log10(mach**5.5 * h * dbar_h * L / r**2)

    spl_blt = g4[cp.newaxis, :, :, :] + g5 + scale[cp.newaxis, :, :, :]

    # Calculation results of p'
    scaled_f_f0 = f_f0 / NSTEP
    scaled_f_f0 = scaled_f_f0[cp.newaxis, :, :, :]

    p6 = cp.sum(10**(spl_blt / 10) * scaled_f_f0, axis=(1, 2, 3))
    p6 /= NPERI

    return p6


def bpm_tip_vs(obj_par, obj_rot, obj_sig, p7):
    # this section should run only for last segment,
    # and the code needs to be fixed.
    NBLADE, NSEG = obj_par.NBLADE, obj_par.NSEG
    NSTEP = obj_par.NSTEP
    NPERI = obj_par.NPERI
    nts = obj_par.nts
    RPM = obj_par.RPM

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    nfreq = obj_par.nfreq

    TE_RCORR = obj_par.TE_RCORR
    C0 = obj_par.C0
    alpha_tip = obj_par.AOA_TIP

    r = obj_rot.r_te
    # L = obj_rot.L
    C = obj_rot.C

    # r[:, :, :-1] = 0.0
    # L[:, :, :-1] = 0.0
    # C[:, :, :-1] = 0.0

    mach = obj_rot.mach
    ufree = obj_rot.ufree

    f_f0 = obj_rot.f_f0_te

    alpha_rat = 2 * cp.pi
    alpha_tipp = alpha_tip * alpha_rat

    # Compute directivity function.
    dbar_h, dbar_l = obj_rot.dbar_h_te, obj_rot.dbar_l_te
    # -----------------------------

    if (TE_RCORR == True):
        L = 0.008 * alpha_tipp * C
    else:
        if (abs(alpha_tipp) <= 2.0):
            L = (0.0230 + 0.0169 * alpha_tipp) * C
        else:
            L = (0.0378 + 0.0095 * alpha_tipp) * C

    mm = (1.0 + 0.036 * alpha_tipp) * mach
    um = mm * C0

    term = mach**2 * mm**3 * L**2 * dbar_h / r**2

    cond = (term != 0.0)
    scale = cp.where(cond, 10.0 * cp.log10(term), 0.0)

    Lum_temp = L / um
    stpp = cen_ob[:, cp.newaxis, cp.newaxis, cp.newaxis
                  ] * Lum_temp[cp.newaxis, :, :, :]

    spl_tip = 26.0 - 30.5 * (cp.log10(stpp) + 0.3)**2 + scale[cp.newaxis,
                                                              :, :, :]

    # Calculation results of p'
    scaled_f_f0 = f_f0 / NSTEP
    scaled_f_f0 = scaled_f_f0[cp.newaxis, :, :, :]

    p7 = cp.sum(10**(spl_tip / 10) * scaled_f_f0, axis=(1, 2, 3))
    p7 /= NPERI

    return p7


def roger_bwi(obj_par, obj_rot, obj_sig, p9):
    NBLADE, NSEG = obj_par.NBLADE, obj_par.NSEG
    NPERI = obj_par.NPERI
    NSTEP = obj_par.NSTEP
    nts = obj_par.nts
    RPM = obj_par.RPM

    # test #
    lti_min = []
    lti_max = []
    # lti = []
    #

    FRQ_BAND = obj_par.FRQ_BAND
    if FRQ_BAND == 'one':
        cen_ob = obj_sig.cen_ob
    elif FRQ_BAND == 'third':
        cen_ob = obj_sig.cen_ob3
    elif FRQ_BAND == 'twelfth':
        cen_ob = obj_sig.cen_ob12

    nfreq = obj_par.nfreq
    C0 = obj_par.C0
    RHO = obj_par.RHO

    re_chord, mach = obj_rot.re_chord, obj_rot.mach
    r, azi, elv = obj_rot.r_le, obj_rot.azi_le, obj_rot.elv_le
    ufree = obj_rot.ufree
    C = obj_rot.C
    L = obj_rot.L
    turb_leng0, turb_ubar = obj_rot.turb_leng, obj_rot.turb_ubar
    # turb_leng0 = obj_rot.turb_leng

    f_f0 = obj_rot.f_f0_le

    # Expand the array's dimension for broadcasting.
    # ----------------------------------------------
    f_f0_ext = f_f0[cp.newaxis, :, :, :]
    ufree_ext = ufree[cp.newaxis, :, :, :]
    C_ext = C[cp.newaxis, :, :, :]
    L_ext = L[cp.newaxis, :, :, :]
    mach_ext = mach[cp.newaxis, :, :, :]
    turb_leng0_ext = turb_leng0[cp.newaxis, cp.newaxis, cp.newaxis, :]
    turb_intensity = turb_ubar[cp.newaxis, cp.newaxis, cp.newaxis, :]

    # re_dh = (1.225 * 0.25 * ufree_ext) / 1.7894e-05
    # turb_intensity = 0.16 * re_dh**(-0.125)
    turb_ubar_ext = turb_intensity * ufree_ext

    # test 12 29
    # mean_vel = cp.mean(ufree, axis=0)
    # fluc_vel = ufree - mean_vel[cp.newaxis, :, :]
    # rms_vel = cp.sqrt(cp.mean(fluc_vel**2, axis=0))
    # # rms shape = (n_blade, n_section)
    # rms_vel_ext = rms_vel[cp.newaxis, cp.newaxis, :, :]
    #

    beta = cp.sqrt(1.0 - mach_ext**2).astype(cp.float64)
    beta2 = beta**2.0

    azir = cp.deg2rad(azi)
    elvr = cp.deg2rad(elv)

    y0 = r * cp.cos(azir)  # along chordwise direction
    x0 = r * cp.sin(azir) * cp.cos(elvr)  # along spanwise direction
    z0 = r * cp.sin(azir) * cp.sin(elvr)  # along rotation axis
    y0_ext = y0[cp.newaxis, :, :, :]
    x0_ext = x0[cp.newaxis, :, :, :]
    z0_ext = z0[cp.newaxis, :, :, :]
    # convection-corrected distance
    s0 = cp.sqrt(y0_ext**2 + beta2 * (x0_ext**2 + z0_ext**2))

    # Compute boundary layer thickness.
    # ---------------------------------
    dstr_p, dstr_s = obj_rot.dstr_p, obj_rot.dstr_s
    dstr_p_ext = dstr_p[cp.newaxis, :, :, :]
    dstr_s_ext = dstr_s[cp.newaxis, :, :, :]

    w_obs = 2.0 * cp.pi * cen_ob
    w_obs_ext = w_obs[:, cp.newaxis, cp.newaxis, cp.newaxis]
    w_src = w_obs_ext / f_f0_ext

    fintg = 0.0
    kac = w_src / C0
    k1 = w_src / ufree_ext

    kstr1 = k1 * (C_ext / 2.0)

    # L_TI computation
    kstr2_min = -(2.0 * L_ext / C_ext) * kac * C_ext / (2.0 * beta)

    kstr2_max = -kstr2_min

    dkstr2 = 0.05 * kac * C_ext / (2.0 * beta)
    dk2 = dkstr2 / (C_ext / 2.0)

    kstr2 = kstr2_min.copy()
    cond = kstr2 + dkstr2 <= kstr2_max

    epsilon = 1e-6

    while cp.any(cond):
        # ===
        k2 = kstr2 / (C_ext / 2.0)

        mu_bar = (kac * C_ext) / (2.0 * beta2)

        kappa_temp = ((kstr2 / beta)**2.0 - mu_bar**2.0).astype(cp.complex128)
        kappa = cp.sqrt(kappa_temp)
        kappai = 1.0j * kappa

        theta0 = kstr1 * mach_ext / (beta * kstr2)
        theta2 = mu_bar * (mach_ext - y0_ext * s0) - cp.pi / 4.0
        theta3 = kappa + 1.0j * mu_bar * (y0_ext / s0)
        # theta4 = 1.0j * theta3
        theta4 = kappai - mu_bar * (y0_ext / s0)

        kappa2 = mu_bar**2.0 * (1.0 / theta0**2.0 - 1.0)
        # turb_length = turb_leng0_ext + dstr_p_ext + dstr_s_ext
        turb_length = turb_leng0_ext
        # kbar = 3/2 * rms_vel_ext**2
        # turb_length = (2 * kbar)**(3/2) / 10**-3

        ke = cp.sqrt(cp.pi) / turb_length *\
            m.gamma(5.0 / 6.0) / m.gamma(1.0 / 3.0)
        khat1 = (k1 / ke).astype(cp.float64)
        khat2 = (k2 / ke).astype(cp.float64)

        # Von karman model
        phiww = 4.0 / (9.0 * cp.pi) * (turb_ubar_ext / ke)**2.0 * \
            (khat1**2.0 + khat2**2.0) / (
                (1.0 + khat1**2.0 + khat2**2.0)**(7.0 / 3.0))

        xarg1 = cp.sqrt((2 * 1.0j * theta3))
        xarg1_cp = abs(xarg1) / cp.sqrt(cp.pi / 2.0)
        # xarg1_cp = xarg1 / cp.sqrt(cp.pi / 2.0)
        xarg1_np = cp.asnumpy(xarg1_cp)
        sfren1, cfren1 = sp.fresnel(xarg1_np)

        diffrac1_np = cfren1 + sfren1 * 1.0j  # Euler's formula
        diffrac1 = cp.asarray(diffrac1_np)

        # Error Function, erf()
        xi_temp = 4.0 * abs(kappa)
        xi = cp.sqrt(xi_temp)

        nt = 101
        nt_temp = cp.arange(nt)
        nt_arr = nt_temp[cp.newaxis, cp.newaxis, cp.newaxis, cp.newaxis, :]
        tt = xi[..., cp.newaxis] / float(nt - 1) * nt_arr
        ferr = cp.sum(cp.exp(-tt[..., 1:-1]**2) * 0.5 *
                      (tt[..., 2:] - tt[..., :-2]), axis=(4))

        ferr *= (2.0 / cp.sqrt(cp.pi))

        # ==================================

        xarg2 = cp.sqrt(2 * (1.0j * kappa + mu_bar * (y0_ext / s0)))
        xarg2_cp = cp.abs(xarg2) / cp.sqrt(cp.pi / 2.0)
        # xarg2_cp = xarg2 / cp.sqrt(cp.pi / 2.0)
        xarg2_np = cp.asnumpy(xarg2_cp)
        sfren2, cfren2 = sp.fresnel(xarg2_np)

        diffrac2_np = cfren2 + 1.0j * sfren2
        diffrac2 = cp.asarray(diffrac2_np)

        xarg3 = cp.sqrt(2.0 * theta4)
        xarg3_cp = cp.abs(xarg3) / cp.sqrt(cp.pi / 2.0)
        # xarg3_cp = xarg3 / cp.sqrt(cp.pi / 2.0)
        xarg3_np = cp.asnumpy(xarg3_cp)
        sfren3, cfren3 = sp.fresnel(xarg3_np)

        diffrac3_np = cfren3 + 1.0j * sfren3
        diffrac3 = cp.asarray(diffrac3_np)

        xarg4 = cp.sqrt(4.0 * kappai)
        xarg4_cp = cp.abs(xarg4) / cp.sqrt(cp.pi / 2.0)
        # xarg4_cp = xarg4 / cp.sqrt(cp.pi / 2.0)
        xarg4_np = cp.asnumpy(xarg4_cp)
        sfren4, cfren4 = sp.fresnel(xarg4_np)

        diffrac4_np = cfren4 + 1.0j * sfren4
        diffrac4 = cp.asarray(diffrac4_np)

        xarg5 = cp.sqrt(2.0 * (kappai + mu_bar * y0_ext / s0))
        xarg5_cp = cp.abs(xarg5) / cp.sqrt(cp.pi / 2.0)
        # xarg5_cp = xarg5 / cp.sqrt(cp.pi / 2.0)
        xarg5_np = cp.asnumpy(xarg5_cp)
        sfren5, cfren5 = sp.fresnel(xarg5_np)

        diffrac5_np = cfren5 + 1.0j * sfren5
        diffrac5 = cp.asarray(diffrac5_np)

        cond_sub = (theta0 < 1.0)  # Sub-critical gust
        cond_sup = ((theta0 > 1.0) & (kappa2 < 0.0))  # Super-critical gust

        num11 = -1.0 / cp.pi
        num12 = cp.sqrt(2.0 / ((kstr1 + 1.0j * beta2 * kappa) * 1.0j * theta3))
        num13 = cp.exp(-1.0j * theta2)
        num14 = cp.sqrt(2.0 / ((kstr1 + beta2 * kappai) * theta4))
        l1 = cp.where(cond_sub,
                      num11 * num12 * num13 * diffrac1,
                      cp.where(cond_sup,
                               num11 * num14 * num13 * diffrac3,
                               0.0))

        num21 = cp.exp(-1.0j * theta2)
        num22 = cp.pi * cp.sqrt(2 * cp.pi * (kstr1 + 1.0j * beta2 * kappa)) * \
            theta3
        num23 = 1.0 - cp.exp(-2.0 * theta3) - ferr
        num24 = 2.0 * cp.exp(-2.0 * theta3) * \
            cp.sqrt(kappa / (1.0j * kappa + mu_bar * (y0_ext / s0)))

        num25 = cp.pi * theta4 * \
            cp.sqrt(2.0 * cp.pi * (kstr1 + beta2 * kappai))
        num26 = 1.0j * (1.0 - cp.exp(2.0j * theta4))
        num27 = cp.exp(2.0j * theta4) * \
            cp.sqrt(2.0 * kappai / (kappai + mu_bar * y0_ext / s0))

        l2 = cp.where(cond_sub,
                      (num21 / num22) * (num23 + num24 * diffrac2),
                      cp.where(cond_sup, (
                          num21 / num25) * (num26 - (1 + 1.0j) *
                                            (diffrac4 - num27 * diffrac5)),
                          0.0))

        # Numerical criteria for L2
        real_Q1 = cp.percentile(l2.real, 25)
        real_Q3 = cp.percentile(l2.real, 75)
        real_IQR = real_Q3 - real_Q1
        real_lower_bound = real_Q1 - 1.5 * real_IQR
        real_upper_bound = real_Q3 + 1.5 * real_IQR

        imag_Q1 = cp.percentile(l2.imag, 25)
        imag_Q3 = cp.percentile(l2.imag, 75)
        imag_IQR = imag_Q3 - imag_Q1
        imag_lower_bound = imag_Q1 - 1.5 * imag_IQR
        imag_upper_bound = imag_Q3 + 1.5 * imag_IQR

        real_outliers = (l2.real < real_lower_bound) \
            | (l2.real > real_upper_bound)
        imag_outliers = (l2.imag < imag_lower_bound) \
            | (l2.imag > imag_upper_bound)

        l2[real_outliers | imag_outliers] = epsilon

        l_ti = l1 + l2

        # l_ti.real[l_ti.real < -100] = -epsilon
        # l_ti.real[l_ti.real > 100] = epsilon
        # l_ti.imag[l_ti.imag < -100] = -epsilon
        # l_ti.imag[l_ti.imag > 100] = epsilon

        # test for output
        lti_min_temp = cp.sign(cp.real(l_ti).min()) * cp.abs(l_ti.min())
        lti_max_temp = cp.sign(cp.real(l_ti).max()) * cp.abs(l_ti.max())
        lti_min_temp[lti_min_temp < -1e+6] = -1e+6
        lti_max_temp[lti_max_temp > 1e+6] = 1e+6
        lti_min.append(lti_min_temp)
        lti_max.append(lti_max_temp)

        fcorr = ((cp.sin((kac * x0_ext / s0 - k2) * (L_ext / 2.0)))**2) / \
            (cp.pi * (L_ext / 2.0) * (kac * x0_ext / s0 - k2)**2.0)

        fintg += (phiww * cp.abs(l_ti)**2.0 * fcorr) * dk2

        # ===
        kstr2 = cp.where(cond, kstr2 + dkstr2, kstr2)
        cond = kstr2 + dkstr2 <= kstr2_max

    spp_psi = ((RHO * kac * C_ext * z0_ext) / (2.0 * s0**2))**2 * cp.pi * \
        ufree_ext * (L_ext / 2.0) * fintg

    # print(f"spp_psi = {spp_psi[:, -1, 0, 0]}")

    # Bandwidth
    bw = cp.zeros_like(cen_ob)

    bw[:-1] = 2.0 * (10.0**((cp.log10(cen_ob[1:]) +
                             cp.log10(cen_ob[:-1])) / 2.0) -
                     cen_ob[:-1])

    bw[-1] = 2.0 * (-10.0**((cp.log10(cen_ob[-1]) +
                             cp.log10(cen_ob[-2])) / 2.0) +
                    cen_ob[-1])

    bw = bw[:, cp.newaxis, cp.newaxis, cp.newaxis]
    cond = (spp_psi == 0.0)
    spl_bwi = cp.where(cond, 1.0e-20,
                       10.0 * cp.log10(spp_psi * bw / (2.0e-5)**2.0))

    # Calculation results of p'
    scaled_f_f0 = f_f0_ext / NSTEP

    p9 = cp.sum(10**(spl_bwi / 10) * scaled_f_f0, axis=(1, 2, 3))
    p9 /= NPERI

    return p9, lti_min, lti_max
    # return p9


# def fresnel(xarg):
#     # Defraction algorithm
#     # Ref.: Michael Möser,
#     #       "Engineering Acoustics: An Introduction to Noise Control,"
#     #       Springer, 2nd Edition, p.342, 2004.
#     x = abs(xarg) / np.sqrt(np.pi / 2.0)
#     # x = xarg
#     arg = np.pi * x**2 / 2.0
#     sin = np.sin(arg)
#     cos = np.cos(arg)

#     cond1 = x > 4.4
#     cond2 = xarg.real < 0.0

#     x4 = x**4
#     x3 = x**3
#     x1 = 0.3183099 - 0.0968 / x4
#     x2 = 0.10132 - 0.154 / x4

#     c_a0 = x.copy()
#     c_somme = x.copy()
#     xmul = -((np.pi / 2.0)**2.0) * (x**4.0)
#     c_an = c_a0.copy()

#     s_a0 = (np.pi / 6.0) * x**3.0
#     s_somme = s_a0.copy()
#     s_an = s_a0.copy()

#     nend = (x + 1.0) * 20.0
#     n = np.zeros_like(nend)
#     cond = n + 1 < nend + 1
#     while np.any(cond):
#         c_xnenn = (2.0 * n + 1.0) * (2.0 * n + 2.0) * (4.0 * n + 5.0)
#         c_an1 = c_an * (4.0 * n + 1.0) * xmul / c_xnenn
#         c_somme += c_an1
#         c_an = c_an1.copy()

#         s_xnenn = (2.0 * n + 2.0) * (2.0 * n + 3.0) * (4.0 * n + 7.0)
#         s_an1 = s_an * (4.0 * n + 3.0) * xmul / s_xnenn
#         s_somme += s_an1
#         s_an = s_an1.copy()

#         n = n + 1
#         cond = n + 1 < nend + 1

#     cfren1 = np.where(cond1, 0.5 + x1 * sin / x - x2 * cos / x3, c_somme)
#     sfren1 = np.where(cond1, 0.5 - x1 * cos / x - x2 * sin / x3, s_somme)
#     cfren1 = np.where(cond2, cfren1 * -1.0, cfren1)
#     sfren1 = np.where(cond2, sfren1 * -1.0, sfren1)

#     return sfren1, cfren1
