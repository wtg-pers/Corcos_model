# import for system.
import os
# import for numerical calculation.
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def apply_rotation_and_get_controlp(obj_par, obj_rot, node):
    NSTEP = obj_par.NSTEP
    NPERI = obj_par.NPERI
    NBLADE = obj_par.NBLADE
    NODEM, NODEN = obj_par.NODEM, obj_par.NODEN

    tk_nseg = obj_rot.tk_nseg
    rot_direc = obj_par.rot_direc

    deform_segs = np.zeros((NSTEP, NBLADE, tk_nseg, 7))

    # Rotate the blade to align with the other blades.
    blades_step = np.radians((360 / NBLADE)) * rot_direc
    for i in range(NBLADE):
        R = np.array([[np.cos(blades_step * i), -np.sin(blades_step * i), 0],
                      [np.sin(blades_step * i), np.cos(blades_step * i), 0],
                      [0, 0, 1]])

        node[:, i, :, :, 0] = (R[0, 0] * node[:, 0, :, :, 0] +
                               R[0, 1] * node[:, 0, :, :, 1] +
                               R[0, 2] * node[:, 0, :, :, 2])
        node[:, i, :, :, 1] = (R[1, 0] * node[:, 0, :, :, 0] +
                               R[1, 1] * node[:, 0, :, :, 1] +
                               R[1, 2] * node[:, 0, :, :, 2])
        node[:, i, :, :, 2] = (R[2, 0] * node[:, 0, :, :, 0] +
                               R[2, 1] * node[:, 0, :, :, 1] +
                               R[2, 2] * node[:, 0, :, :, 2])

    # Rotate the blades to align with the time-series position.
    psi_step = int(360 / NSTEP)
    psi = [i for i in range(0, 360, psi_step)]
    psi = np.array(np.deg2rad(psi)) * rot_direc

    ii = 0
    for it in psi:
        Rz = np.array([[np.cos(it), -np.sin(it), 0],
                      [np.sin(it), np.cos(it), 0],
                      [0, 0, 1]])

        node[ii, :, :, :, 0] = (Rz[0, 0] * node[0, :, :, :, 0] +
                                Rz[0, 1] * node[0, :, :, :, 1] +
                                Rz[0, 2] * node[0, :, :, :, 2])
        node[ii, :, :, :, 1] = (Rz[1, 0] * node[0, :, :, :, 0] +
                                Rz[1, 1] * node[0, :, :, :, 1] +
                                Rz[1, 2] * node[0, :, :, :, 2])
        node[ii, :, :, :, 2] = (Rz[2, 0] * node[0, :, :, :, 0] +
                                Rz[2, 1] * node[0, :, :, :, 1] +
                                Rz[2, 2] * node[0, :, :, :, 2])

        ii += 1

    # Calculate the control(center)-points & normal vectors.
    controlp = 0.25 * (node[:, :, :NODEN-1, :NODEM-1, :] +
                       node[:, :, :NODEN-1, 1:NODEM, :] +
                       node[:, :, 1:NODEN, 1:NODEM, :] +
                       node[:, :, 1:NODEN, :NODEM-1, :])

    vector1 = node[:, :, 1:NODEN, 1:NODEM, :] - \
        node[:, :, :NODEN-1, :NODEM-1, :]

    vector2 = node[:, :, 1:NODEN, :NODEM-1, :] - \
        node[:, :, :NODEN-1, 1:NODEM, :]

    normal_x = (vector1[..., 1] * vector2[..., 2] -
                vector1[..., 2] * vector2[..., 1])
    normal_y = (-vector1[..., 0] * vector2[..., 2] +
                vector1[..., 2] * vector2[..., 0])
    normal_z = (vector1[..., 0] * vector2[..., 1] -
                vector1[..., 1] * vector2[..., 0])

    normal_vector = np.array([normal_x, normal_y, normal_z])

    # Normalize the normal vector
    abs_v = np.sqrt(np.einsum('ijklm,ijklm->jklm',
                              normal_vector, normal_vector))
    normal_vector /= abs_v

    # # Adjust normal vector direction to always point outward
    # z_axis = np.array([0, 0, 1])
    # ref_points = np.array([0, 0, 0])
    # ref_points_ext = ref_points[np.newaxis, np.newaxis,
    #                             np.newaxis, np.newaxis, :]
    # ref_vector = controlp[..., :] - ref_points_ext

    # direc_vector = np.einsum('jklmi, i -> ijklm',
    #                          ref_vector, z_axis)

    # # dot_prod = np.einsum('ijklm, jklmi -> ijklm', normal_vector, ref_vector)
    # # normal_vector = np.where(dot_prod < 0, -normal_vector, normal_vector)
    # normal_vector = np.where(direc_vector < 0,
    #                          -normal_vector, normal_vector)

    # Data array re-shaping to transfer in ojbect.
    controlp = np.reshape(controlp, (NSTEP, NBLADE, tk_nseg, 3))

    normal_vector = np.transpose(normal_vector, (1, 2, 3, 4, 0))
    normal_vector = np.reshape(normal_vector, (NSTEP, NBLADE, tk_nseg, 3))

    abs_v = np.reshape(abs_v, (NSTEP, NBLADE, tk_nseg))

    # Data transfer.
    deform_segs[:, :, :, :3] = controlp
    deform_segs[:, :, :, 3:6] = normal_vector
    deform_segs[:, :, :, 6] = 0.5 * abs_v

    deform_segs = np.tile(deform_segs, (NPERI, 1, 1, 1))

    obj_rot.deform_segs = cp.asarray(deform_segs)

    nv_temp = os.path.join(obj_par.dir_out, "normal_vector.dat")
    fnv = open(nv_temp, 'w')
    fnv.write("""variables="x", "y", "z", "nx", "ny", "nz"\n""")
    for k in range(deform_segs.shape[1]):
        fnv.write(f'zone t="blade{k+1}", i={NODEM-1}, j={NODEN-1}\n')
        for j in range(deform_segs.shape[2]):
            fnv.write(
                f"{deform_segs[0, k, j, 0]}\t"
                f"{deform_segs[0, k, j, 1]}\t"
                f"{deform_segs[0, k, j, 2]}\t"
                f"{deform_segs[0, k, j, 3]}\t"
                f"{deform_segs[0, k, j, 4]}\t"
                f"{deform_segs[0, k, j, 5]}\n")
    fnv.close()

    return None


def spherical_metrix(obj_par, obj_rot, temp_mics):
    NSTEP, NBLADE = obj_par.NSTEP, obj_par.NBLADE
    rot_direc = obj_par.rot_direc

    temp_mics = cp.asnumpy(temp_mics)

    x_le, y_le, z_le = obj_rot.x_le, obj_rot.y_le, obj_rot.z_le
    x_te, y_te, z_te = obj_rot.x_te, obj_rot.y_te, obj_rot.z_te

    alpha_geor = cp.asnumpy(cp.deg2rad(obj_rot.alpha_geo))

    x_le = cp.asnumpy(x_le)
    y_le = cp.asnumpy(y_le)
    z_le = cp.asnumpy(z_le)

    x_te = cp.asnumpy(x_te)
    y_te = cp.asnumpy(y_te)
    z_te = cp.asnumpy(z_te)

    # De-rotate L.E & T.E coordinates.
    # =================================
    psi_step = int(360 / NSTEP)
    psi = [i for i in range(0, 360, psi_step)]
    psi = np.array(np.deg2rad(psi)) * rot_direc
    mrot = -1 * rot_direc

    ii = 0
    for it in psi:
        Rz = np.array([[np.cos(it * mrot), -np.sin(it * mrot), 0],
                      [np.sin(it * mrot), np.cos(it * mrot), 0],
                      [0, 0, 1]])

        x_le[ii, :, :] = (Rz[0, 0] * x_le[0, :, :] +
                          Rz[0, 1] * y_le[0, :, :] +
                          Rz[0, 2] * z_le[0, :, :])
        y_le[ii, :, :] = (Rz[1, 0] * x_le[0, :, :] +
                          Rz[1, 1] * y_le[0, :, :] +
                          Rz[1, 2] * z_le[0, :, :])
        z_le[ii, :, :] = (Rz[2, 0] * x_le[0, :, :] +
                          Rz[2, 1] * y_le[0, :, :] +
                          Rz[2, 2] * z_le[0, :, :])

        x_te[ii, :, :] = (Rz[0, 0] * x_te[0, :, :] +
                          Rz[0, 1] * y_te[0, :, :] +
                          Rz[0, 2] * z_te[0, :, :])
        y_te[ii, :, :] = (Rz[1, 0] * x_te[0, :, :] +
                          Rz[1, 1] * y_te[0, :, :] +
                          Rz[1, 2] * z_te[0, :, :])
        z_te[ii, :, :] = (Rz[2, 0] * x_te[0, :, :] +
                          Rz[2, 1] * y_te[0, :, :] +
                          Rz[2, 2] * z_te[0, :, :])

        ii += 1
    # =================================

    #### 09-11 by JSH ####
    # De-rotate along the angle of blade number
    # Rotate coordinates according at Omega.
    # ======================================
    blade_step = np.radians((360 / NBLADE)) * rot_direc
    for i in range(NBLADE):
        Rz = np.array([[np.cos(-blade_step * i), -np.sin(-blade_step * i), 0],
                      [np.sin(-blade_step * i), np.cos(-blade_step * i), 0],
                      [0, 0, 1]])

        x_le[:, i, :] = (Rz[0, 0] * x_le[:, 0, :] +
                         Rz[0, 1] * y_le[:, 0, :] +
                         Rz[0, 2] * z_le[:, 0, :])
        y_le[:, i, :] = (Rz[1, 0] * x_le[:, 0, :] +
                         Rz[1, 1] * y_le[:, 0, :] +
                         Rz[1, 2] * z_le[:, 0, :])
        z_le[:, i, :] = (Rz[2, 0] * x_le[:, 0, :] +
                         Rz[2, 1] * y_le[:, 0, :] +
                         Rz[2, 2] * z_le[:, 0, :])

        x_te[:, i, :] = (Rz[0, 0] * x_te[:, 0, :] +
                         Rz[0, 1] * y_te[:, 0, :] +
                         Rz[0, 2] * z_te[:, 0, :])
        y_te[:, i, :] = (Rz[1, 0] * x_te[:, 0, :] +
                         Rz[1, 1] * y_te[:, 0, :] +
                         Rz[1, 2] * z_te[:, 0, :])
        z_te[:, i, :] = (Rz[2, 0] * x_te[:, 0, :] +
                         Rz[2, 1] * y_te[:, 0, :] +
                         Rz[2, 2] * z_te[:, 0, :])
    # ======================================
    ####

    # Transform coordinates for centering segment to be zero.
    # =======================================================
    rx_le = temp_mics[0] - x_le
    ry_le = temp_mics[1] - y_le
    rz_le = temp_mics[2] - z_le
    r_le = np.sqrt(rx_le**2 + ry_le**2 + rz_le**2)

    rx_te = temp_mics[0] - x_te
    ry_te = temp_mics[1] - y_te
    rz_te = temp_mics[2] - z_te
    r_te = np.sqrt(rx_te**2 + ry_te**2 + rz_te**2)

    # =======================================================

    # Rotate L.E & T.E coordinates according at Omega.
    # =================================
    psi_step = int(360 / NSTEP)
    psi = [i for i in range(0, 360, psi_step)]
    psi = np.array(np.deg2rad(psi)) * rot_direc

    ii = 0
    for it in psi:
        Rz = np.array([[np.cos(it), -np.sin(it), 0],
                      [np.sin(it), np.cos(it), 0],
                      [0, 0, 1]])

        rx_le[ii, :, :] = (Rz[0, 0] * rx_le[0, :, :] +
                           Rz[0, 1] * ry_le[0, :, :] +
                           Rz[0, 2] * rz_le[0, :, :])
        ry_le[ii, :, :] = (Rz[1, 0] * rx_le[0, :, :] +
                           Rz[1, 1] * ry_le[0, :, :] +
                           Rz[1, 2] * rz_le[0, :, :])
        rz_le[ii, :, :] = (Rz[2, 0] * rx_le[0, :, :] +
                           Rz[2, 1] * ry_le[0, :, :] +
                           Rz[2, 2] * rz_le[0, :, :])

        rx_te[ii, :, :] = (Rz[0, 0] * rx_te[0, :, :] +
                           Rz[0, 1] * ry_te[0, :, :] +
                           Rz[0, 2] * rz_te[0, :, :])
        ry_te[ii, :, :] = (Rz[1, 0] * rx_te[0, :, :] +
                           Rz[1, 1] * ry_te[0, :, :] +
                           Rz[1, 2] * rz_te[0, :, :])
        rz_te[ii, :, :] = (Rz[2, 0] * rx_te[0, :, :] +
                           Rz[2, 1] * ry_te[0, :, :] +
                           Rz[2, 2] * rz_te[0, :, :])

        ii += 1
    # =================================
    # ======================================
    blades_step = np.radians((360 / NBLADE)) * rot_direc
    for i in range(NBLADE):
        Rz = np.array([[np.cos(blades_step * i), -np.sin(blades_step * i), 0],
                      [np.sin(blades_step * i), np.cos(blades_step * i), 0],
                      [0, 0, 1]])

        rx_le[:, i, :] = (Rz[0, 0] * rx_le[:, 0, :] +
                          Rz[0, 1] * ry_le[:, 0, :] +
                          Rz[0, 2] * rz_le[:, 0, :])
        ry_le[:, i, :] = (Rz[1, 0] * rx_le[:, 0, :] +
                          Rz[1, 1] * ry_le[:, 0, :] +
                          Rz[1, 2] * rz_le[:, 0, :])
        rz_le[:, i, :] = (Rz[2, 0] * rx_le[:, 0, :] +
                          Rz[2, 1] * ry_le[:, 0, :] +
                          Rz[2, 2] * rz_le[:, 0, :])

        rx_te[:, i, :] = (Rz[0, 0] * rx_te[:, 0, :] +
                          Rz[0, 1] * ry_te[:, 0, :] +
                          Rz[0, 2] * rz_te[:, 0, :])
        ry_te[:, i, :] = (Rz[1, 0] * rx_te[:, 0, :] +
                          Rz[1, 1] * ry_te[:, 0, :] +
                          Rz[1, 2] * rz_te[:, 0, :])
        rz_te[:, i, :] = (Rz[2, 0] * rx_te[:, 0, :] +
                          Rz[2, 1] * ry_te[:, 0, :] +
                          Rz[2, 2] * rz_te[:, 0, :])
    # ======================================

    # Rotate coordinates according to geometric angle.
    # ================================================
    temp1_le = rx_le
    temp2_le = ry_le
    temp3_le = rz_le

    rx_le = temp1_le
    ry_le = temp2_le * np.cos(-alpha_geor) - temp3_le * np.sin(-alpha_geor)
    rz_le = temp2_le * np.sin(-alpha_geor) + temp3_le * np.cos(-alpha_geor)

    temp1_te = rx_te
    temp2_te = ry_te
    temp3_te = rz_te

    rx_te = temp1_te
    ry_te = temp2_te * np.cos(-alpha_geor) - temp3_te * np.sin(-alpha_geor)
    rz_te = temp2_te * np.sin(-alpha_geor) + temp3_te * np.cos(-alpha_geor)

    hxz_le = np.hypot(rx_le, rz_le)
    hxz_te = np.hypot(rx_te, rz_te)

    # ==== 09.27 by JSH ====
    elv_le = np.rad2deg(np.arctan2(rz_le, rx_le))
    azi_le = np.rad2deg(np.arctan2(hxz_le, -ry_le))

    elv_te = np.rad2deg(np.arctan2(rz_te, rx_te))
    azi_te = np.rad2deg(np.arctan2(hxz_te, -ry_te))

    # elv = np.rad2deg(np.arctan2(rz, rx))
    # azi = np.rad2deg(np.arctan2(np.sqrt(rx*rx + rz*rz), -ry))

    # te_x = cp.asarray(te_x)
    # te_y = cp.asarray(te_y)
    # te_z = cp.asarray(te_z)
    # rx = cp.asarray(rx)
    # ry = cp.asarray(ry)
    # rz = cp.asarray(rz)

    obj_rot.r_le = cp.asarray(r_le)
    obj_rot.elv_le = cp.asarray(elv_le)
    obj_rot.azi_le = cp.asarray(azi_le)

    obj_rot.r_te = cp.asarray(r_te)
    obj_rot.elv_te = cp.asarray(elv_te)
    obj_rot.azi_te = cp.asarray(azi_te)

    return None
