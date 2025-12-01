function G = Gravity ( m , Vg , dtij , dRij , rcm )

    G = m * Vg' * (dtij + dRij * rcm);

end