
    th(1) = q(1,mm);
    th(2) = q(2,mm);
    r(3) = q(3,mm);
    th(4) = q(4,mm);
    qp = dq(:,mm);
    qpp = ddq(:,mm);

    for k = 1 : n
        T(:,:,k) = T_Matrix(q(k),alph(k),d(k),r(k));
    end

    Ti(:,:,1) = T(:,:,1);
    for k = 2 : n
        Ti(:,:,k) = Ti(:,:,k-1) * T(:,:,k);
    end

    zi(:,1:n) = Ti(1:3,3,1:n);
    ti(:,1:n) = Ti(1:3,4,1:n);
    Ri(:,:,1:n) = Ti(1:3,1:3,1:n);
 
    for i = 1 : n
        for k = 1 : i
            li(:,k,i) = ti(:,i) - ti(:,k);
            [ Jvi(:,k,i) , Jwi(:,k,i) ] = Jacobian_Matrix ( zi(:,k) , li(:,k,i) );
        end
    end

    for i = 1 : n
        Jgi(:,:,i) = Jgi_Matrix ( Jvi(:,:,i) , Ri(:,:,i) , rcm_skew(:,:,i) , Jwi(:,:,i) );
    end
    Jgi = simplify(Jgi); %%%%

    M = Inertia_Matrix ( m(1) , Jgi(:,:,1) , Jwi(:,:,1) , Ri(:,:,1) , Ig(:,:,1) );
    for k = 2 : n
        M = M + Inertia_Matrix ( m(k) , Jgi(:,:,k) , Jwi(:,:,k) , Ri(:,:,k) , Ig(:,:,k) );
    end
    M = simplify(M);    %%%%

    for k = 1 : n
        dT(:,:,k) = dT_Matrix ( q(k) , alph(k) , type(k) );
    end

    for i = 1 : n
        for j = 1 : n
            if j > i
                dTij(:,:,i,j) = zeros(4);
            elseif j == i
                if j == 1
                    dTij(:,:,i,j) = dT(:,:,j);
                else
                    dTij(:,:,i,j) = Ti(:,:,j-1) * dT(:,:,j);
                end
            else
                if j == 1 
                    Tij = eye(4);
                    for k = j + 1 : i
                        Tij = Tij * Ti(:,:,k);
                    end
                    dTij(:,:,i,j) = dT(:,:,j) * Tij;
                else
                    Tij = eye(4);
                    for k = j + 1 : i
                        Tij = Tij * Ti(:,:,k);
                    end
                    dTij(:,:,i,j) = Ti(:,:,j-1) * dT(:,:,j) * Tij;
                end
            end
            dRij(:,:,i,j) = dTij(1:3,1:3,i,j);
            dzij(:,i,j) = dTij(1:3,3,i,j);
            dtij(:,i,j) = dTij(1:3,4,i,j);
        end
    end

    for j = 1 : n
        G(j,1) = Gravity ( m(1) , Vg , dtij(:,1,j) , dRij(:,:,1,j) , rcm(:,1) );
        for i = 2 : n
            G(j,1) = G(j,1) + Gravity ( m(i) , Vg , dtij(:,i,j) , dRij(:,:,i,j) , rcm(:,i) );
        end
    end

    for i = 1 : n
        for k = 1 : i
            for j = 1 : n
                dlij(:,k,i,j) = dtij(:,i,j) - dtij(:,k,j);
                [ dJvi(:,k,i,j) , dJwi(:,k,i,j) ] = dJ_Matrix ( zi(:,k) , li(:,k,i) , dzij(:,k,j) , dlij(:,k,i,j) );
            end
        end
    end

    for i = 1 : n
        for j = 1 : n
            dJgi(:,:,i,j) = dJgi_Matrix ( dJvi(:,:,i,j) , dRij(:,:,i,j) , rcm_skew(:,:,i) , Ri(:,:,i) , Jwi(:,:,i) , dJwi(:,:,i,j) );
        end
    end

    for j = 1 : n
        dMj(:,:,j) = dMj_Matrix ( m(1) , dJgi(:,:,1,j) , Jgi(:,:,1) , dJwi(:,:,1,j) , Ri(:,:,1) , Ig(:,:,1) , Jwi(:,:,1) , dRij(:,:,1,j) );
        for i = 2 : n
            dMj(:,:,j) = dMj(:,:,j) + dMj_Matrix ( m(i) , dJgi(:,:,i,j) , Jgi(:,:,i) , dJwi(:,:,i,j) , Ri(:,:,i) , Ig(:,:,i) , Jwi(:,:,i) , dRij(:,:,i,j) );
        end
    end

    Mp = Mp_Matrix ( dMj(:,:,1) , qp(1) );
    for j = 2 : n
        Mp = Mp + Mp_Matrix ( dMj(:,:,j) , qp(j) );
    end

    for j = 1 : n
        N2(j,:) = qp' * dMj(:,:,j);
    end
  
    N = Mp - 1/2 * N2;

    Tau_id = simplify(vpa(M * qpp + N * qp + G - Fv .* qp - Fc .* sign(qp)));


