function dJgi = dJgi_Matrix ( dJvi , dRij , rcm_skew , Ri , Jwi , dJwi )

    dJgi = dJvi - dRij * rcm_skew * (Ri)' * Jwi - Ri * rcm_skew * (dRij)' * Jwi - Ri * rcm_skew * (Ri)' * dJwi;

end