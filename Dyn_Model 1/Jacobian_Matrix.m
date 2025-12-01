function [ Jvi , Jwi ] = Jacobian_Matrix ( zi , li )

    Jvi = cross(zi,li);
    Jwi = zi;

end