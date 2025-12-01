function dMj = dMj_Matrix ( m , dJgi , Jgi , dJwi , Ri , Ig , Jwi , dRij )

    dMj = m * ( (dJgi)' * Jgi + (Jgi)' * dJgi ) + (dJwi)' * Ri * Ig * (Ri)' * Jwi + (Jwi)' * dRij * Ig * (Ri)' * Jwi ...
    + (Jwi)' * Ri * Ig * (dRij)' * Jwi + (Jwi)' * Ri * Ig * (Ri)' * dJwi;

end