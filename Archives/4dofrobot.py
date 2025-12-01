#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fonctionnel donne le meme resultat qu'en matlab mais de façon plus longue 
import time
import sympy as sp
import numpy as np
from scipy.io import loadmat, savemat
from multiprocessing import Pool, cpu_count
# ---------------------------------------------------------
#  Fonction skew : matrice antisymétrique associée à un vecteur
# ---------------------------------------------------------
def skew(v):
    """
    Retourne la matrice antisymétrique 3x3 associée au vecteur 3x1 v.
    v : sympy.Matrix (3, 1)
    """
    return sp.Matrix([
        [0,      -v[2],  v[1]],
        [v[2],   0,     -v[0]],
        [-v[1],  v[0],   0]
    ])


# ---------------------------------------------------------
#  Script principal (traduction de dynamic_model.m)
# ---------------------------------------------------------
if __name__ == "__main__":

    start_time = time.time()

    # -----------------------------
    #  Paramètres généraux
    # -----------------------------
    n = 4
    joint_type = [1, 1, 0, 1]   # type = [1;1;0;1];

    # -----------------------------
    #  Paramètres MDH du TS2-40
    # -----------------------------
    d    = sp.Matrix([0, 0.22, 0.24, 0])
    alph = sp.Matrix([0, 0,    0,    0])
    r    = sp.Matrix([0, 0,    0,    0])
    th   = sp.Matrix([0, 0,    0,    0])

    # Gravité
    Vg = sp.Matrix([0, 0, 9.81])

    # -----------------------------
    #  Paramètres dynamiques symboliques
    # -----------------------------
    m  = sp.Matrix(sp.symbols("m1:%d"  % (n+1),  real=True))

    xx = sp.Matrix(sp.symbols("xx1:%d" % (n+1), real=True))
    yy = sp.Matrix(sp.symbols("yy1:%d" % (n+1), real=True))
    zz = sp.Matrix(sp.symbols("zz1:%d" % (n+1), real=True))
    xy = sp.Matrix(sp.symbols("xy1:%d" % (n+1), real=True))
    xz = sp.Matrix(sp.symbols("xz1:%d" % (n+1), real=True))
    yz = sp.Matrix(sp.symbols("yz1:%d" % (n+1), real=True))

    Ig = []
    for k in range(n):
        Ig_k = sp.Matrix([
            [xx[k], xy[k], xz[k]],
            [xy[k], yy[k], yz[k]],
            [xz[k], yz[k], zz[k]]
        ])
        Ig.append(Ig_k)

    x = sp.Matrix(sp.symbols("x1:%d" % (n+1), real=True))
    y = sp.Matrix(sp.symbols("y1:%d" % (n+1), real=True))
    z = sp.Matrix(sp.symbols("z1:%d" % (n+1), real=True))

    rcm = sp.zeros(3, n)
    for k in range(n):
        rcm[:, k] = sp.Matrix([x[k], y[k], z[k]])

    rcm_skew = [skew(rcm[:, k]) for k in range(n)]

    Fc = sp.Matrix(sp.symbols("Fc1:%d" % (n+1), real=True))
    Fv = sp.Matrix(sp.symbols("Fv1:%d" % (n+1), real=True))

    # ---------------------------------------------------------
    #  Chargement des données .mat (qf, dqf, ddqf, tauf, q)
    # ---------------------------------------------------------
    print("Chargement des fichiers .mat ...")

    qf_full   = loadmat("qf.mat")["qf"]       # 6 x N
    dqf_full  = loadmat("dqf.mat")["dqf"]     # 6 x N
    ddqf_full = loadmat("ddqf.mat")["ddqf"]   # 6 x N
    tauf_full = loadmat("tauf.mat")["tauf"]   # 6 x N

    # On garde seulement les 4 premiers axes
    qf   = qf_full[0:n, :]
    dqf  = dqf_full[0:n, :]
    ddqf = ddqf_full[0:n, :]
    tauf = tauf_full[0:n, :]

    q_data = loadmat("q.mat")
    t = np.squeeze(q_data["t"])  # 1 x N

    # ---------------------------------------------------------
    #  Sous-échantillonnage (un échantillon sur 380)
    # ---------------------------------------------------------
    print("Sous-échantillonnage des données ...")

    nn = int((len(t) - 1) / 380)

    t_sub   = np.zeros((1, nn))
    q_sub   = np.zeros((n, nn))
    dq_sub  = np.zeros((n, nn))
    ddq_sub = np.zeros((n, nn))
    tau_sub = np.zeros((n, nn))

    for k in range(nn):
        idx = k * 380
        t_sub[0, k]   = t[idx]
        q_sub[:, k]   = qf[:, idx]
        dq_sub[:, k]  = dqf[:, idx]
        ddq_sub[:, k] = ddqf[:, idx]
        tau_sub[:, k] = tauf[:, idx]

    # ---------------------------------------------------------
    #  Préallocation de tau_id (symbolique)
    # ---------------------------------------------------------
    tau_id = sp.MutableDenseMatrix.zeros(n, nn)

    print(f"Nombre d'échantillons après sous-échantillonnage : nn = {nn}")

    # Pour ne pas exploser le temps de calcul, on peut limiter ici
    # comme dans le commentaire MATLAB :
    #  % for mm = 1:10 si pas réussi a faire l'optimisation
    MAX_SAMPLES =   10001        # mets 10 ici si tu veux tester plus vite
    mm_max = min(nn, MAX_SAMPLES)

    print(f"Début du calcul dynamique pour {mm_max} échantillons ...")

    # ---------------------------------------------------------
    #  Boucle sur les échantillons
    # ---------------------------------------------------------
    for mm in range(mm_max):

        pct = 100.0 * (mm + 1) / mm_max
        print(f"\n===== Échantillon {mm+1}/{mm_max} ({pct:5.1f} %) =====")

        # ---------------------------------------------
        #  Données pour l'échantillon mm
        # ---------------------------------------------
        th[0] = float(q_sub[0, mm])
        th[1] = float(q_sub[1, mm])
        r[2]  = float(q_sub[2, mm])
        th[3] = float(q_sub[3, mm])

        qp  = sp.Matrix(dq_sub[:, mm])
        qpp = sp.Matrix(ddq_sub[:, mm])

        # ---------------------------------------------
        #  Matrices de transformation T(:,:,k)
        # ---------------------------------------------
        T = []
        for k in range(n):
            ck = sp.cos(th[k])
            sk = sp.sin(th[k])
            ca = sp.cos(alph[k])
            sa = sp.sin(alph[k])

            Tk = sp.Matrix([
                [ck,        -sk,         0,         d[k]],
                [ca * sk,   ca * ck,    -sa, -r[k] * sa],
                [sa * sk,   sa * ck,     ca,  r[k] * ca],
                [0,          0,          0,         1]
            ])
            T.append(Tk)

        # Ti(:,:,k)
        Ti = [None] * n
        Ti[0] = T[0]
        for k in range(1, n):
            Ti[k] = Ti[k-1] * T[k]

        # zi, ti, Ri
        zi = []
        ti_vec = []
        Ri = []
        for k in range(n):
            Ti_k = Ti[k]
            Ri_k = Ti_k[0:3, 0:3]
            ti_k = Ti_k[0:3, 3]
            zi_k = Ti_k[0:3, 2]
            Ri.append(Ri_k)
            ti_vec.append(ti_k)
            zi.append(zi_k)

        # ---------------------------------------------
        #  Jacobiennes Jvi, Jwi
        # ---------------------------------------------
        Jvi = []
        Jwi = []
        li = {}  # li[(i,k)]

        for i in range(n):
            Jvi_i = sp.zeros(3, n)
            Jwi_i = sp.zeros(3, n)
            for k in range(i+1):
                li_ik = ti_vec[i] - ti_vec[k]
                li[(i, k)] = li_ik
                Jvi_i[:, k] = zi[k].cross(li_ik)
                Jwi_i[:, k] = zi[k]
            Jvi.append(Jvi_i)
            Jwi.append(Jwi_i)

        # Jgi(:,:,i)
        Jgi = []
        for i in range(n):
            Jgi_i = Jvi[i] - Ri[i] * rcm_skew[i] * Ri[i].T * Jwi[i]
            Jgi.append(Jgi_i)

        # ---------------------------------------------
        #  Matrice d'inertie M
        # ---------------------------------------------
        M = (m[0] * (Jgi[0].T * Jgi[0]) +
             Jwi[0].T * Ri[0] * Ig[0] * Ri[0].T * Jwi[0])

        for k in range(1, n):
            M += (m[k] * (Jgi[k].T * Jgi[k]) +
                  Jwi[k].T * Ri[k] * Ig[k] * Ri[k].T * Jwi[k])

        # ---------------------------------------------
        #  Dérivée des matrices de transformation
        # ---------------------------------------------
        dT = []
        for k in range(n):
            ck = sp.cos(th[k])
            sk = sp.sin(th[k])
            ca = sp.cos(alph[k])
            sa = sp.sin(alph[k])

            if joint_type[k] == 1:  # rotatif
                dTk = sp.Matrix([
                    [-sk,          -ck,        0, 0],
                    [ca * ck, -ca * sk,        0, 0],
                    [sa * ck, -sa * sk,        0, 0],
                    [0,            0,          0, 1]
                ])
            else:  # prismatique
                dTk = sp.Matrix([
                    [0, 0, 0, 0],
                    [0, 0, 0, -sa],
                    [0, 0, 0,  ca],
                    [0, 0, 0, 1]
                ])
            dT.append(dTk)

        # dTij, dRij, dzij, dtij
        dTij = {}
        dRij = {}
        dzij = {}
        dtij = {}

        for i in range(n):
            for j in range(n):
                if j > i:
                    dTij[(i, j)] = sp.zeros(4)
                elif j == i:
                    if j == 0:
                        dTij[(i, j)] = dT[j]
                    else:
                        dTij[(i, j)] = Ti[j-1] * dT[j]
                else:  # j < i
                    Tij_tmp = sp.eye(4)
                    for kk in range(j+1, i+1):
                        Tij_tmp = Tij_tmp * Ti[kk]
                    if j == 0:
                        dTij[(i, j)] = dT[j] * Tij_tmp
                    else:
                        dTij[(i, j)] = Ti[j-1] * dT[j] * Tij_tmp

                dRij[(i, j)] = dTij[(i, j)][0:3, 0:3]
                dzij[(i, j)] = dTij[(i, j)][0:3, 2]
                dtij[(i, j)] = dTij[(i, j)][0:3, 3]

        # ---------------------------------------------
        #  Vecteur de gravité G
        # ---------------------------------------------
        G = sp.Matrix.zeros(n, 1)
        for j in range(n):
            term = m[0] * (Vg.T *
                           (dtij[(0, j)] + dRij[(0, j)] * rcm[:, 0]))[0, 0]
            for i in range(1, n):
                term += m[i] * (Vg.T *
                                (dtij[(i, j)] + dRij[(i, j)] * rcm[:, i]))[0, 0]
            G[j, 0] = term

        # ---------------------------------------------
        #  Dérivées des Jacobiennes
        # ---------------------------------------------
        dlij = {}
        dJvi_dict = {}
        dJwi_dict = {}

        for i in range(n):
            for k in range(i+1):
                for j in range(n):
                    dl = dtij[(i, j)] - dtij[(k, j)]
                    dlij[(i, k, j)] = dl
                    dJvi_dict[(i, k, j)] = (
                        dzij[(k, j)].cross(li[(i, k)]) +
                        zi[k].cross(dl)
                    )
                    dJwi_dict[(i, k, j)] = dzij[(k, j)]

        dJgi = {}
        dJvi_full = {}
        dJwi_full = {}

        for i in range(n):
            for j in range(n):
                dJvi_ij = sp.zeros(3, n)
                dJwi_ij = sp.zeros(3, n)
                for k in range(i+1):
                    dJvi_ij[:, k] = dJvi_dict[(i, k, j)]
                    dJwi_ij[:, k] = dJwi_dict[(i, k, j)]

                dJvi_full[(i, j)] = dJvi_ij
                dJwi_full[(i, j)] = dJwi_ij

                dJgi[(i, j)] = (
                    dJvi_ij
                    - dRij[(i, j)] * rcm_skew[i] * Ri[i].T * Jwi[i]
                    - Ri[i] * rcm_skew[i] * dRij[(i, j)].T * Jwi[i]
                    - Ri[i] * rcm_skew[i] * Ri[i].T * dJwi_ij
                )

        # ---------------------------------------------
        #  dMj(:,:,j)
        # ---------------------------------------------
        dMj = {}

        for j in range(n):
            dMj_j = (
                m[0] *
                (dJgi[(0, j)].T * Jgi[0] + Jgi[0].T * dJgi[(0, j)]) +
                dJwi_full[(0, j)].T * Ri[0] * Ig[0] * Ri[0].T * Jwi[0] +
                Jwi[0].T * dRij[(0, j)] * Ig[0] * Ri[0].T * Jwi[0] +
                Jwi[0].T * Ri[0] * Ig[0] * dRij[(0, j)].T * Jwi[0] +
                Jwi[0].T * Ri[0] * Ig[0] * Ri[0].T * dJwi_full[(0, j)]
            )

            for i in range(1, n):
                dMj_j += (
                    m[i] *
                    (dJgi[(i, j)].T * Jgi[i] + Jgi[i].T * dJgi[(i, j)]) +
                    dJwi_full[(i, j)].T * Ri[i] * Ig[i] * Ri[i].T * Jwi[i] +
                    Jwi[i].T * dRij[(i, j)] * Ig[i] * Ri[i].T * Jwi[i] +
                    Jwi[i].T * Ri[i] * Ig[i] * dRij[(i, j)].T * Jwi[i] +
                    Jwi[i].T * Ri[i] * Ig[i] * Ri[i].T * dJwi_full[(i, j)]
                )

            dMj[j] = dMj_j

        # Mp
        Mp = dMj[0] * qp[0]
        for j in range(1, n):
            Mp += dMj[j] * qp[j]

        # N2 et N
        N2 = sp.Matrix.zeros(n, n)
        for j in range(n):
            N2[j, :] = (qp.T * dMj[j])

        N = Mp - sp.Rational(1, 2) * N2

        # ---------------------------------------------
        #  Couple tau_id
        # ---------------------------------------------
                # ---------------------------------------------
        #  Couple tau_id (version simplifiée)
        # ---------------------------------------------
        visc    = Fv.multiply_elementwise(qp)
        coulomb = Fc.multiply_elementwise(qp.applyfunc(sp.sign))

        tau_mm = M * qpp + N * qp + G - visc - coulomb  # expression symbolique (avec Floats)

        # Liste des variables “intéressantes” pour le collect
        vars_poly = list(m) + list(x) + list(y) + list(z) + list(zz) + list(Fc) + list(Fv)

        tau_mm_simpl = sp.Matrix.zeros(n, 1)
        for i in range(n):
            expr = tau_mm[i]

            # 1) Développer
            expr = sp.expand(expr)

            # 2) Regrouper par les variables dynamiques / géométriques
            expr = sp.collect(expr, vars_poly)

            # 3) Optionnel : légère simplification supplémentaire
            expr = sp.simplify(expr)

            # 4) Optionnel : arrondir proprement les coefficients
            expr = sp.nsimplify(expr, [*vars_poly], rational=True)

            tau_mm_simpl[i] = expr

        # On stocke la version simplifiée (symbolique) dans tau_id
        tau_id[:, mm] = tau_mm_simpl

        elapsed = time.time() - start_time
        print(f"  -> échantillon {mm+1} terminé, temps total ~ {elapsed:6.1f} s")


    print("\nSauvegarde de tau_sub.mat et tau_id.mat ...")

    # tau_sub : numérique, OK
    savemat("tau_sub.mat", {"tau_sub": tau_sub})

    tau_id_str = np.empty((n, mm_max), dtype=object)
    for i in range(n):
        for j in range(mm_max):
            tau_id_str[i, j] = str(tau_id[i, j])
    savemat("tau_idpython.mat", {"tau_idpython": tau_id_str})
    total_time = time.time() - start_time
    print(f"Terminé en {total_time:.1f} secondes pour {mm_max} échantillons.")

