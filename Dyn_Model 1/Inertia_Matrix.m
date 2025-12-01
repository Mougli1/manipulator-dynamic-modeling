function M =  Inertia_Matrix ( m , Jgi , Jwi , Ri , Ig )

    M = m * (Jgi)' * Jgi + (Jwi)' * Ri * Ig * (Ri)' * Jwi;

end