function Jgi = Jgi_Matrix ( Jvi , Ri , rcm_skew , Jwi )

Jgi = Jvi - Ri * rcm_skew * (Ri)' * Jwi;

end