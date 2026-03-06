# import for system
import os
import sys
import time
from tqdm import tqdm
import numpy as np
import cupy as cp
from numba import cuda
import matplotlib.pyplot as plt


def lowson(obj_par, obj_rot, temp_mics):
    NBLADE = obj_par.NBLADE
    C0 = obj_par.C0
    mach_fwd = obj_par.mach_fwd
    omega = obj_par.omega
    nts, dt = obj_par.nts, obj_par.dt

    blades = [0]

    nseg = obj_par.NODEN - 1
    p_new = obj_rot.ls_p_new

    j_start, j_end = obj_par.INTERP_LVL, (obj_par.nts - 2)

    p = cp.zeros((nts, NBLADE, nseg))

    # ====
    threads_per_blk = (32, 1, 16)
    blk_x = int(cp.ceil(p_new.shape[0] / threads_per_blk[0]))
    blk_y = int(cp.ceil((p_new.shape[1]) / threads_per_blk[1]))
    blk_z = int(cp.ceil((p_new.shape[2]) / threads_per_blk[2]))
    blk_per_grid = (blk_x, blk_y, blk_z)
    # ====

    lowson_segs = obj_rot.lowson_segs
    x0 = lowson_segs[:, :, :, 0]
    y0 = lowson_segs[:, :, :, 1]
    z0 = lowson_segs[:, :, :, 2]
    fx = lowson_segs[:, :, :, 3]
    fy = lowson_segs[:, :, :, 4]
    fz = lowson_segs[:, :, :, 5]

    dx0 = temp_mics[0] - x0[j_start-1:j_end, :, :]
    dy0 = temp_mics[1] - y0[j_start-1:j_end, :, :]
    dz0 = temp_mics[2] - z0[j_start-1:j_end, :, :]

    diffs = cp.stack((dx0, dy0, dz0), axis=-1)

    R0 = cp.linalg.norm(diffs, axis=-1)

    # Generate a source time for contribution of tonal pressure.
    src_t = cp.einsum('IJK, I -> IJK',
                      cp.ones((nts, NBLADE, nseg)),
                      dt * cp.arange(nts))

    rx = -1 * dx0
    obs_ang = cp.arccos(rx / R0)

    wind_eff = (mach_fwd * cp.cos(obs_ang) +
                cp.sqrt(1 - mach_fwd**2 * cp.sin(obs_ang)**2)
                ) / (1 - mach_fwd**2)

    obs_t = cp.zeros((nts, NBLADE, nseg))

    obs_t[j_start-1:j_end, :, :] = (src_t[j_start-1:j_end, :, :] +
                                    R0 / C0 * wind_eff)

    obj_rot.ld_src_t = src_t
    obj_rot.ld_obs_t = obs_t

    aa = cp.sqrt(x0[j_start-1:j_end, :, :]**2 +
                 y0[j_start-1:j_end, :, :]**2)

    cosgamma = x0[j_start-1:j_end, :, :] / aa
    singamma = y0[j_start-1:j_end, :, :] / aa

    M = aa * omega / C0
    Mx = -M * singamma - mach_fwd
    My = M * cosgamma
    Mz = cp.zeros_like(M)

    dMx = -M * cosgamma * omega
    dMy = -M * singamma * omega
    dMz = cp.zeros_like(dMx)

    r = cp.sqrt(dx0**2 + dy0**2 + dz0**2)
    rx = dx0 / r
    ry = dy0 / r
    rz = dz0 / r

    Mr = Mx * rx + My * ry + Mz * rz
    dMr = dMx * rx + dMy * ry + dMz * rz

    fx_tmp = fx[j_start-2:j_end+1, :, :]
    fy_tmp = fy[j_start-2:j_end+1, :, :]
    fz_tmp = fz[j_start-2:j_end+1, :, :]

    dlx = (fx_tmp[2:, :, :] - fx_tmp[:-2, :, :]) / (2 * dt)
    dly = (fy_tmp[2:, :, :] - fy_tmp[:-2, :, :]) / (2 * dt)
    dlz = (fz_tmp[2:, :, :] - fz_tmp[:-2, :, :]) / (2 * dt)

    fx_cent = fx_tmp[1:-1, :, :]
    fy_cent = fy_tmp[1:-1, :, :]
    fz_cent = fz_tmp[1:-1, :, :]
    lr = fx_cent * rx + fy_cent * ry + fz_cent * rz
    dlr = dlx * rx + dly * ry + dlz * rz
    liMi = fx_cent * Mx + fy_cent * My + fz_cent * Mz

    coeff1 = 1/4/cp.pi

    p_near = (coeff1 * (1.0 / r**2 / (1.0 - Mr)) *
              (lr * (1.0 - M**2) / (1.0 - Mr) - liMi))
    p_far = (coeff1 * (1.0 / r / C0 / (1.0 - Mr)
                       ) * (dlr + lr * dMr / (1.0 - Mr)))

    p[j_start-1:j_end, :, :] = p_near + p_far

    # get new_t
    get_new_t(obj_par, obj_rot, "loading")

    obj_rot.ld_p = p
    new_t = obj_rot.ld_new_t

    # conduct interpolation
    lagrange_poly_interp[blk_per_grid, threads_per_blk](new_t, obs_t,
                                                        p, p_new,
                                                        j_start, j_end)

    # print(p_new)
    obj_rot.ld_p_new = p_new

    sum_sound_pressure(obj_par, obj_rot, "loading")

    return None


def function_farassat_1a(obj_par, obj_rot, temp_mics, key):

    nmics = obj_par.nmics
    OM, ON = obj_par.OM, obj_par.ON
    NBLADE = obj_par.NBLADE
    RHO = obj_par.RHO
    C0 = obj_par.C0
    mach_fwd = obj_par.mach_fwd
    omega = obj_par.omega
    nts, dt = obj_par.nts, obj_par.dt

    blades = [0]

    if key == "loading":
        nseg = obj_rot.ld_nseg
        p_new = obj_rot.ld_p_new
    elif key == "thickness":
        nseg = obj_rot.tk_nseg
        p_new = obj_rot.tk_p_new

    j_start, j_end = obj_par.INTERP_LVL, (obj_par.nts - 2)

    p = cp.zeros((nts, NBLADE, nseg))

    # temp_p_new = cp.zeros((nts, NBLADE, ld_nseg))

    # ====
    threads_per_blk = (32, 1, 16)
    blk_x = int(cp.ceil(p_new.shape[0] / threads_per_blk[0]))
    blk_y = int(cp.ceil((p_new.shape[1]) / threads_per_blk[1]))
    blk_z = int(cp.ceil((p_new.shape[2]) / threads_per_blk[2]))
    blk_per_grid = (blk_x, blk_y, blk_z)
    # ====

    if key == "loading":
        cen_segs = obj_rot.cen_segs
        x0 = cen_segs[:, :, :, 0]
        y0 = cen_segs[:, :, :, 1]
        z0 = cen_segs[:, :, :, 2]
        fx = cen_segs[:, :, :, 3]
        fy = cen_segs[:, :, :, 4]
        fz = cen_segs[:, :, :, 5]
        area = cen_segs[:, :, :, 6]
    elif key == "thickness":
        deform_segs = obj_rot.deform_segs
        x0 = deform_segs[:, :, :, 0]
        y0 = deform_segs[:, :, :, 1]
        z0 = deform_segs[:, :, :, 2]
        ux = deform_segs[:, :, :, 3]
        uy = deform_segs[:, :, :, 4]
        uz = deform_segs[:, :, :, 5]
        area = deform_segs[:, :, :, 6]

    dx0 = temp_mics[0] - x0[j_start-1:j_end, :, :]
    dy0 = temp_mics[1] - y0[j_start-1:j_end, :, :]
    dz0 = temp_mics[2] - z0[j_start-1:j_end, :, :]

    diffs = cp.stack((dx0, dy0, dz0), axis=-1)

    R0 = cp.linalg.norm(diffs, axis=-1)

    # Generate a source time for contribution of tonal pressure.
    src_t = cp.einsum('IJK, I -> IJK',
                      cp.ones((nts, NBLADE, nseg)),
                      dt * cp.arange(nts))

    # for wind effect
    rx = -1 * dx0
    obs_ang = cp.arccos(rx / R0)

    wind_eff = (mach_fwd * cp.cos(obs_ang) +
                cp.sqrt(1 - mach_fwd**2 * cp.sin(obs_ang)**2)
                ) / (1 - mach_fwd**2)

    obs_t = cp.zeros((nts, NBLADE, nseg))

    obs_t[j_start-1:j_end, :, :] = (src_t[j_start-1:j_end, :, :] +
                                    R0 / C0 * wind_eff)

    if key == "loading":

        obj_rot.ld_src_t = src_t
        obj_rot.ld_obs_t = obs_t

        aa = cp.sqrt(x0[j_start-1:j_end, :, :]**2 +
                     y0[j_start-1:j_end, :, :]**2)

        cosgamma = x0[j_start-1:j_end, :, :] / aa
        singamma = y0[j_start-1:j_end, :, :] / aa

        M = aa * omega / C0
        Mx = -M * singamma - mach_fwd
        My = M * cosgamma
        Mz = cp.zeros_like(M)

        dMx = -M * cosgamma * omega
        dMy = -M * singamma * omega
        dMz = cp.zeros_like(dMx)

        r = cp.sqrt(dx0**2 + dy0**2 + dz0**2)
        rx = dx0 / r
        ry = dy0 / r
        rz = dz0 / r

        Mr = Mx * rx + My * ry + Mz * rz
        dMr = dMx * rx + dMy * ry + dMz * rz

        fx_tmp = fx[j_start-2:j_end+1, :, :]
        fy_tmp = fy[j_start-2:j_end+1, :, :]
        fz_tmp = fz[j_start-2:j_end+1, :, :]

        dlx = (fx_tmp[2:, :, :] - fx_tmp[:-2, :, :]) / (2 * dt)
        dly = (fy_tmp[2:, :, :] - fy_tmp[:-2, :, :]) / (2 * dt)
        dlz = (fz_tmp[2:, :, :] - fz_tmp[:-2, :, :]) / (2 * dt)

        fx_cent = fx_tmp[1:-1, :, :]
        fy_cent = fy_tmp[1:-1, :, :]
        fz_cent = fz_tmp[1:-1, :, :]
        lr = fx_cent * rx + fy_cent * ry + fz_cent * rz
        dlr = dlx * rx + dly * ry + dlz * rz
        liMi = fx_cent * Mx + fy_cent * My + fz_cent * Mz

        coeff1 = 1/4/cp.pi

        p_near = (coeff1 / C0 * (dlr / r / (1-Mr)**2) *
                  area[j_start-1:j_end, :, :] + coeff1 *
                  ((lr - liMi) / r**2 / (1-Mr)**2) *
                  area[j_start-1:j_end, :, :])

        p_far = (coeff1 / C0 *
                 (lr * (r * dMr + C0 * Mr - C0 * M**2) / r**2 /
                  (1-Mr)**3) * area[j_start-1:j_end, :, :])

        p[j_start-1:j_end, :, :] = p_near + p_far

        # get new_t
        get_new_t(obj_par, obj_rot, "loading")

        obj_rot.ld_p = p
        new_t = obj_rot.ld_new_t

        # conduct interpolation
        lagrange_poly_interp[blk_per_grid, threads_per_blk](new_t, obs_t,
                                                            p, p_new,
                                                            j_start, j_end)

        obj_rot.ld_p_new = p_new

        sum_sound_pressure(obj_par, obj_rot, "loading")

    elif key == "thickness":

        obj_rot.tk_src_t = src_t
        obj_rot.tk_obs_t = obs_t

        amx = x0[j_start-2:j_end-1, :, :]
        ax = x0[j_start-1:j_end, :, :]
        apx = x0[j_start:j_end+1, :, :]
        Mx = (apx - amx) / (2.0 * dt) / C0 - mach_fwd
        dMx = (apx - 2.0 * ax + amx) / dt**2 / C0

        amy = y0[j_start-2:j_end-1, :, :]
        ay = y0[j_start-1:j_end, :, :]
        apy = y0[j_start:j_end+1, :, :]
        My = (apy - amy) / (2.0 * dt) / C0
        dMy = (apy - 2.0 * ay + amy) / dt**2 / C0

        amz = z0[j_start-2:j_end-1, :, :]
        az = z0[j_start-1:j_end, :, :]
        apz = z0[j_start:j_end+1, :, :]
        Mz = (apz - amz) / (2.0 * dt) / C0
        dMz = (apz - 2.0 * az + amz) / dt**2 / C0

        M = cp.sqrt(Mx**2 + My**2 + Mz**2)
        r = cp.sqrt(dx0**2 + dy0**2 + dz0**2)
        rx = dx0 / r
        ry = dy0 / r
        rz = dz0 / r

        Mr = Mx * rx + My * ry + Mz * rz
        dMr = dMx * rx + dMy * ry + dMz * rz

        ux_tmp = ux[j_start-2:j_end+1, :, :]
        uy_tmp = uy[j_start-2:j_end+1, :, :]
        uz_tmp = uz[j_start-2:j_end+1, :, :]

        Mn = Mx * ux_tmp[1:-1, :, :] + My * uy_tmp[1:-1, :, :] + \
            Mz * uz_tmp[1:-1, :, :]
        dMn = dMx * ux_tmp[1:-1, :, :] + dMy * uy_tmp[1:-1, :, :] + \
            dMz * uz_tmp[1:-1, :, :]

        dux = (ux_tmp[2:, :, :] - ux_tmp[:-2, :, :]) / (2*dt)
        duy = (uy_tmp[2:, :, :] - uy_tmp[:-2, :, :]) / (2*dt)
        duz = (uz_tmp[2:, :, :] - uz_tmp[:-2, :, :]) / (2*dt)

        Mdn = Mx * dux + My * duy + Mz * duz

        coeff1 = 1 / 4 / cp.pi * RHO

        p_thick1 = (coeff1 * (dMn * C0 + Mdn * C0) / r / (1 - Mr)**2 *
                    area[j_start-1:j_end, :, :])
        p_thick2 = (coeff1 * C0 * Mn *
                    (r * dMr + C0 * Mr - C0 * M**2) /
                    r**2 / (1 - Mr)**3 * area[j_start-1:j_end, :, :])
        p[j_start-1:j_end, :, :] = p_thick1 + p_thick2

        # get new_t
        get_new_t(obj_par, obj_rot, "thickness")

        obj_rot.tk_p = p
        new_t = obj_rot.tk_new_t

        # conduct interpolation
        lagrange_poly_interp[blk_per_grid, threads_per_blk](new_t, obs_t,
                                                            p, p_new,
                                                            j_start, j_end)

        obj_rot.tk_p_new = p_new

        sum_sound_pressure(obj_par, obj_rot, "thickness")

    return None


def get_new_t(obj_par, obj_rot, key):
    j_start, j_end = obj_par.INTERP_LVL, (obj_par.nts - 2)
    NINTERP = obj_par.NINTERP

    if key == "loading":
        obs_t, new_t = obj_rot.ld_obs_t, obj_rot.ld_new_t
    elif key == "thickness":
        obs_t, new_t = obj_rot.tk_obs_t, obj_rot.tk_new_t

    new_t[0] = cp.amax(obs_t[j_start - 1, :, :])
    new_t[NINTERP - 1] = cp.amin(obs_t[j_end - 1, :, :])

    dt_ip = (new_t[NINTERP - 1] - new_t[0]) / (NINTERP - 1)

    for i in range(1, NINTERP):
        new_t[i] = new_t[0] + dt_ip * float(i)

    obj_rot.dt_ip = dt_ip

    if key == "loading":
        obj_rot.ld_new_t = new_t
    elif key == "thickness":
        obj_rot.tk_new_t = new_t

    # print(f"dt_ip shape = {dt_ip.shape}")
    # print(f"new_t shape = {new_t.shape}")

    return None


@cuda.jit
def lagrange_poly_interp(new_t, t_ob, p, p_new, j_start, j_end):
    i, j, k = cuda.grid(3)
    if i < p_new.shape[0] and j < p_new.shape[1] and k < p_new.shape[2]:

        x_in = new_t[i]

        for jj in range(j_start - 1, j_end):

            if x_in >= t_ob[jj, j, k] and x_in <= t_ob[jj + 1, j, k]:
                if jj == 1:
                    x_0, x_1, x_2, x_3 = t_ob[jj:jj+4, j, k]
                    f_l0, f_l1, f_l2, f_l3 = p[jj:jj+4, j, k]
                elif jj == j_end - 1:
                    x_0, x_1, x_2, x_3 = t_ob[jj-3:jj+1, j, k]
                    f_l0, f_l1, f_l2, f_l3 = p[jj-3:jj+1, j, k]
                else:
                    x_0, x_1, x_2, x_3 = t_ob[jj-1:jj+3, j, k]
                    f_l0, f_l1, f_l2, f_l3 = p[jj-1:jj+3, j, k]

                # Lagrange polynomial interpolation
                den = (x_0 - x_1) * (x_0 - x_2) * (x_0 - x_3)
                in_x0 = (x_in - x_1) * (x_in - x_2) * (x_in - x_3) / den

                den = (x_1 - x_0) * (x_1 - x_2) * (x_1 - x_3)
                in_x1 = (x_in - x_0) * (x_in - x_2) * (x_in - x_3) / den

                den = (x_2 - x_0) * (x_2 - x_1) * (x_2 - x_3)
                in_x2 = (x_in - x_0) * (x_in - x_1) * (x_in - x_3) / den

                den = (x_3 - x_0) * (x_3 - x_1) * (x_3 - x_2)
                in_x3 = (x_in - x_0) * (x_in - x_1) * (x_in - x_2) / den

                p_new[i, j, k] = in_x0 * f_l0 + in_x1 * \
                    f_l1 + in_x2 * f_l2 + in_x3 * f_l3

                break


def sum_sound_pressure(obj_par, obj_rot, key):
    NINTERP = obj_par.NINTERP

    if key == "loading":
        p_new = obj_rot.ld_p_new
    elif key == "thickness":
        p_new = obj_rot.tk_p_new

    p_sum = cp.sum(p_new, axis=(1, 2))

    aver_p = cp.sum(p_sum) / NINTERP

    temp = p_sum - aver_p

    if key == "loading":
        # obj_rot.ld_p_sum = cp.asarray(temp)
        obj_rot.ld_p_sum = temp
    elif key == "thickness":
        # obj_rot.tk_p_sum = cp.asarray(temp)
        obj_rot.tk_p_sum = temp

    return None


@cuda.jit
def synthesis_tonal_noise(ld_time, tk_time, new_t2,
                          NINTERP, ld_p, ld_p2, tk_p, tk_p2):
    i, j = cuda.grid(2)

    if i == 0 and j == 0:

        new_t2[0] = max(ld_time[0], tk_time[0])

        new_t2[-1] = min(ld_time[-1], tk_time[-1])

        dt_ip2 = (new_t2[-1] - new_t2[0]) / (NINTERP - 1)
        for idx in range(1, NINTERP):
            new_t2[idx] = new_t2[0] + dt_ip2 * float(idx)

    cuda.syncthreads()

    if i < NINTERP - 1 and j < NINTERP:
        for idx in range(NINTERP - 1):
            if (new_t2[j] >= ld_time[idx] and new_t2[j] <= ld_time[idx + 1]):
                ld_p2[j] = ((ld_p[idx + 1] - ld_p[idx]) /
                            (ld_time[idx + 1] - ld_time[idx]) *
                            (new_t2[j] - ld_time[idx]) +
                            ld_p[idx])
                break

        for idx in range(NINTERP - 1):
            if (new_t2[j] >= tk_time[idx] and new_t2[j] <= tk_time[idx + 1]):
                tk_p2[j] = ((tk_p[idx + 1] - tk_p[idx]) /
                            (tk_time[idx + 1] - tk_time[idx]) *
                            (new_t2[j] - tk_time[idx]) +
                            tk_p[idx])
                break


def do_synthesis(obj_par, obj_rot, new_t2, ld_p2, tk_p2, tonal_p):
    NINTERP = obj_par.NINTERP
    ld_time = obj_rot.ld_new_t
    tk_time = obj_rot.tk_new_t
    ld_p_sum = obj_rot.ld_p_sum
    tk_p_sum = obj_rot.tk_p_sum
    # ====
    threads_per_blk = (16, 16)
    blk_x = int(cp.ceil(NINTERP / threads_per_blk[0]))
    blk_y = int(cp.ceil(NINTERP / threads_per_blk[1]))
    blk_per_grid = (blk_x, blk_y)
    # ====

    synthesis_tonal_noise[blk_per_grid, threads_per_blk](ld_time, tk_time,
                                                         new_t2, NINTERP,
                                                         ld_p_sum, ld_p2,
                                                         tk_p_sum, tk_p2)

    temp2 = ld_p2[:] + tk_p2[:]

    aver_p = 0.0
    aver_p = sum(temp2) / NINTERP

    tonal_p[:] = temp2[:] - aver_p
    ld_p2[:] = ld_p2[:] - aver_p
    tk_p2[:] = tk_p2[:] - aver_p

    return None


def doppler_comp(obj_par, obj_rot, temp_mics):
    # Input data's coordinate have a information about time-step, blade N.
    NSTEP = obj_par.NSTEP
    NPERI = obj_par.NPERI
    NBLADE = obj_par.NBLADE
    NSEG = obj_par.NSEG
    nts = obj_par.nts
    rot_direc = obj_par.rot_direc
    RPM = obj_par.RPM
    C0 = obj_par.C0

    x_le, y_le, z_le = obj_rot.x_le, obj_rot.y_le, obj_rot.z_le
    x_te, y_te, z_te = obj_rot.x_te, obj_rot.y_te, obj_rot.z_te
    a_star, a_geo = obj_rot.alpha_star, obj_rot.alpha_geo

    a_ind = a_star - a_geo  # induced angle of attack

    temp_mics = cp.asnumpy(temp_mics)

    x_le = cp.asnumpy(x_le)
    y_le = cp.asnumpy(y_le)
    z_le = cp.asnumpy(z_le)
    x_te = cp.asnumpy(x_te)
    y_te = cp.asnumpy(y_te)
    z_te = cp.asnumpy(z_te)
    alpha_ind = cp.asnumpy(a_ind)
    alpha_indr = np.deg2rad(alpha_ind)

    rx_le = temp_mics[0] - x_le
    ry_le = temp_mics[1] - y_le
    rz_le = temp_mics[2] - z_le
    dr_le = np.sqrt(rx_le**2 + ry_le**2 + rz_le**2)
    r_le = np.sqrt(x_le**2 + y_le**2 + z_le**2)

    rx_te = temp_mics[0] - x_te
    ry_te = temp_mics[1] - y_te
    rz_te = temp_mics[2] - z_te
    dr_te = np.sqrt(rx_te**2 + ry_te**2 + rz_te**2)
    r_te = np.sqrt(x_te**2 + y_te**2 + z_te**2)

    evx = np.zeros((nts+1, NBLADE, NSEG), dtype=cp.float64)
    evy = np.ones((nts+1, NBLADE, NSEG), dtype=cp.float64)
    evz = np.zeros((nts+1, NBLADE, NSEG), dtype=cp.float64)

    # Rotate ev vector along induced angle of attack.
    # ===============================================
    temp1 = evx.copy()
    temp2 = evy.copy()
    temp3 = evz.copy()

    evx = temp1
    evy = temp2 * np.cos(alpha_indr) - temp3 * np.sin(alpha_indr)
    evz = temp2 * np.sin(alpha_indr) + temp3 * np.cos(alpha_indr)
    # ===============================================

    # Rotate ev vector along azimuth angle.
    # =====================================
    temp1 = evx
    temp2 = evy

    blade_step = np.radians((360 / NBLADE)) * rot_direc
    for i in range(NBLADE):
        Rz = np.array([[np.cos(-blade_step * i), np.sin(-blade_step * i), 0],
                      [-np.sin(-blade_step * i), np.cos(-blade_step * i), 0],
                      [0, 0, 1]])
        if i != 0:
            evx[:, i, :] = (Rz[0, 0] * temp1[:, 0, :] +
                            Rz[0, 1] * temp2[:, 0, :])

            evy[:, i, :] = (Rz[1, 0] * temp1[:, 0, :] +
                            Rz[1, 1] * temp2[:, 0, :])

    psi_step = int(360 / NSTEP)
    psi = [i for i in range(0, (360 * NPERI), psi_step)]
    psi = np.array(np.deg2rad(psi)) * rot_direc

    ii = 0
    for it in psi:
        Rz = np.array([[np.cos(-it), np.sin(-it), 0],
                      [-np.sin(-it), np.cos(-it), 0],
                      [0, 0, 1]])

        if it != 0:
            evx[ii, :, :] = (Rz[0, 0] * evx[0, :, :] +
                             Rz[0, 1] * evy[0, :, :])

            evy[ii, :, :] = (Rz[1, 0] * evx[0, :, :] +
                             Rz[1, 1] * evy[0, :, :])

        ii += 1

    evx[-1, :, :] = evx[0, :, :].copy()
    evy[-1, :, :] = evy[0, :, :].copy()
    evz[-1, :, :] = evz[0, :, :].copy()

    # ======================================

    r_dot_er_le = rx_le * evx + ry_le * evy + rz_le * evz
    cos_ksi_le = r_dot_er_le / dr_le

    r_dot_er_te = rx_te * evx + ry_te * evy + rz_te * evz
    cos_ksi_te = r_dot_er_te / dr_te

    omega_dot = (RPM / 60) * 2 * np.pi

    f_f0_le = 1.0 - omega_dot * r_le * cos_ksi_le / C0

    f_f0_te = 1.0 - omega_dot * r_te * cos_ksi_te / C0

    cos_ksi_le = cp.asarray(cos_ksi_le)
    f_f0_le = cp.asarray(f_f0_le)

    cos_ksi_te = cp.asarray(cos_ksi_te)
    f_f0_te = cp.asarray(f_f0_te)

    obj_rot.cos_ksi_le = cos_ksi_le
    obj_rot.f_f0_le = f_f0_le

    obj_rot.cos_ksi_te = cos_ksi_te
    obj_rot.f_f0_te = f_f0_te

    return None


def calculation_spl(p):
    cond = (p != 0.0)
    spl = cp.where(cond, 10.0 * cp.log10(p), 0.0)

    return spl
