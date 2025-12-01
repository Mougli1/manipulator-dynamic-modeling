function [ dJvi , dJwi ] = dJ_Matrix ( zi , li , dzij , dlij )

    dJvi = cross(dzij,li) + cross(zi,dlij);
    dJwi = dzij;

end