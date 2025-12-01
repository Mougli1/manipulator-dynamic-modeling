clear
close all
clc

%% === 1. INITIALISATION ===
n = 4;
d = [ 0 ; 0.22 ; 0.24 ; 0 ];
alph = [ 0 ; 0 ; 0 ; 0 ];
r = [ 0 ; 0 ; 0 ; 0 ];
th = [ 0 ; 0 ; 0 ; 0 ];
Vg = [ 0 ; 0 ; 9.81 ];
digits(4);
type = [ 1 ; 1 ; 0 ; 1 ];

m = sym("m",[ n 1 ],'real');
xx = sym("xx",[ n 1 ],'real');
yy = sym("yy",[ n 1 ],'real');
zz = sym("zz",[ n 1 ],'real');
xy = sym("xy",[ n 1 ],'real');
xz = sym("xz",[ n 1 ],'real');
yz = sym("yz",[ n 1 ],'real');
for k = 1 : n
    Ig(:,:,k) = [ xx(k) xy(k) xz(k) ; xy(k) yy(k) yz(k) ; xz(k) yz(k) zz(k) ];
end

x = sym("x",[ n 1 ],'real');
y = sym("y",[ n 1 ],'real');
z = sym("z",[ n 1 ],'real');
rcm = [ x' ; y' ; z' ];
for k = 1 : n
    rcm_skew(:,:,k) = skew (rcm(:,k));
end
Fc = sym("Fc",[ n 1 ],'real');
Fv = sym("Fv",[ n 1 ],'real');

%% === 2. CHARGEMENT ===
load q.mat
load qf.mat 
load dqf.mat 
load ddqf.mat
load tauf.mat

nn = floor((length(t)-1)/380); 

% Reconstruction des variables comme dans le Main du prof
for k = 1 : nn
    tt(1,k) = t(1,(k-1)*350+1);
    q(:,k) = qf(:,(k-1)*350+1);
    dq(:,k) = dqf(:,(k-1)*350+1);
    ddq(:,k) = ddqf(:,(k-1)*350+1);
    tau_sub_sampled(:,k) = tauf(:,(k-1)*350+1);
end

Tau_id = sym(zeros(n,nn)); % Initialisation (qui sera inutile à cause de l'écrasement)

tic
fprintf('Calcul en cours...\n');

%% === 3. BOUCLE (Comportement strict du prof) ===
for mm = 1 : 4
    if mod(mm, 10) == 0, disp(mm); end

    % --- CONTENU DE LOOP.M ---
    th(1) = q(1,mm);
    th(2) = q(2,mm);
    r(3) = q(3,mm);
    th(4) = q(4,mm);
    qp = dq(:,mm);
    qpp = ddq(:,mm);

    for k = 1 : n
        % Utilisation de q(k) comme le prof (fige l'angle à t=0 pour k=1,2,4, erreur potentielle mais fidèle au code)
        T(:,:,k) = T_Matrix_Internal(q(k),alph(k),d(k),r(k));
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
            [ Jvi(:,k,i) , Jwi(:,k,i) ] = Jacobian_Matrix_Internal ( zi(:,k) , li(:,k,i) );
        end
    end

    for i = 1 : n
        Jgi(:,:,i) = Jgi_Matrix_Internal ( Jvi(:,:,i) , Ri(:,:,i) , rcm_skew(:,:,i) , Jwi(:,:,i) );
    end
    Jgi = simplify(Jgi); 

    M = Inertia_Matrix_Internal ( m(1) , Jgi(:,:,1) , Jwi(:,:,1) , Ri(:,:,1) , Ig(:,:,1) );
    for k = 2 : n
        M = M + Inertia_Matrix_Internal ( m(k) , Jgi(:,:,k) , Jwi(:,:,k) , Ri(:,:,k) , Ig(:,:,k) );
    end
    M = simplify(M);    

    for k = 1 : n
        dT(:,:,k) = dT_Matrix_Internal ( q(k) , alph(k) , type(k) );
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
        G(j,1) = Gravity_Internal ( m(1) , Vg , dtij(:,1,j) , dRij(:,:,1,j) , rcm(:,1) );
        for i = 2 : n
            G(j,1) = G(j,1) + Gravity_Internal ( m(i) , Vg , dtij(:,i,j) , dRij(:,:,i,j) , rcm(:,i) );
        end
    end

    for i = 1 : n
        for k = 1 : i
            for j = 1 : n
                dlij(:,k,i,j) = dtij(:,i,j) - dtij(:,k,j);
                [ dJvi(:,k,i,j) , dJwi(:,k,i,j) ] = dJ_Matrix_Internal ( zi(:,k) , li(:,k,i) , dzij(:,k,j) , dlij(:,k,i,j) );
            end
        end
    end

    for i = 1 : n
        for j = 1 : n
            dJgi(:,:,i,j) = dJgi_Matrix_Internal ( dJvi(:,:,i,j) , dRij(:,:,i,j) , rcm_skew(:,:,i) , Ri(:,:,i) , Jwi(:,:,i) , dJwi(:,:,i,j) );
        end
    end

    for j = 1 : n
        dMj(:,:,j) = dMj_Matrix_Internal ( m(1) , dJgi(:,:,1,j) , Jgi(:,:,1) , dJwi(:,:,1,j) , Ri(:,:,1) , Ig(:,:,1) , Jwi(:,:,1) , dRij(:,:,1,j) );
        for i = 2 : n
            dMj(:,:,j) = dMj(:,:,j) + dMj_Matrix_Internal ( m(i) , dJgi(:,:,i,j) , Jgi(:,:,i) , dJwi(:,:,i,j) , Ri(:,:,i) , Ig(:,:,i) , Jwi(:,:,i) , dRij(:,:,i,j) );
        end
    end

    Mp = Mp_Matrix_Internal ( dMj(:,:,1) , qp(1) );
    for j = 2 : n
        Mp = Mp + Mp_Matrix_Internal ( dMj(:,:,j) , qp(j) );
    end

    for j = 1 : n
        N2(j,:) = qp' * dMj(:,:,j);
    end
  
    N = Mp - 1/2 * N2;

    % === MODIFICATION ICI POUR COLLER AU PROF ===
    % PAS d'index (:,mm) => On écrase la variable Tau_id à chaque tour.
    % À la fin, Tau_id sera un vecteur 4x1 (valeur de la dernière itération).
    Tau_id = simplify(vpa(M * qpp + N * qp + G - Fv .* qp - Fc .* sign(qp)));

end
save Tau_id.mat Tau_id
save tau_sub_sampled.mat tau_sub_sampled
toc


%% === FONCTIONS LOCALES ===

function T = T_Matrix_Internal(th,alph,d,r)
    T = [ cos(th)             , - sin(th)           , 0          , d                  ; 
          cos(alph) * sin(th) , cos(alph) * cos(th) , -sin(alph) , - r * sin(alph) ;
          sin(alph) * sin(th) , sin(alph) * cos(th) , cos(alph)  , r * cos(alph)   ;             
          0                   , 0                   , 0          , 1                 ];
end

function [ Jvi , Jwi ] = Jacobian_Matrix_Internal ( zi , li )
    Jvi = cross(zi,li);
    Jwi = zi;
end

function Jgi = Jgi_Matrix_Internal ( Jvi , Ri , rcm_skew , Jwi )
    Jgi = Jvi - Ri * rcm_skew * (Ri)' * Jwi;
end

function M = Inertia_Matrix_Internal ( m , Jgi , Jwi , Ri , Ig )
    M = m * (Jgi)' * Jgi + (Jwi)' * Ri * Ig * (Ri)' * Jwi;
end

function dT = dT_Matrix_Internal ( th , alph , type )
    if type == 1
        dT = [ - sin(th)             , - cos(th)             , 0 , 0 ; 
                 cos(alph) * cos(th) , - cos(alph) * sin(th) , 0 , 0 ;
                 sin(alph) * cos(th) , - sin(alph) * sin(th) , 0 , 0 ;             
                 0                   , 0                     , 0 , 1 ];
    elseif type == 0
        dT = [ 0 , 0 , 0 , 0 ; 
               0 , 0 , 0 , - sin(alph) ;
               0 , 0 , 0 , cos(alph) ;             
               0 , 0 , 0 , 1 ];
    end
end

function G = Gravity_Internal ( m , Vg , dtij , dRij , rcm )
    G = m * Vg' * (dtij + dRij * rcm);
end

function [ dJvi , dJwi ] = dJ_Matrix_Internal ( zi , li , dzij , dlij )
    dJvi = cross(dzij,li) + cross(zi,dlij);
    dJwi = dzij;
end

function dJgi = dJgi_Matrix_Internal ( dJvi , dRij , rcm_skew , Ri , Jwi , dJwi )
    dJgi = dJvi - dRij * rcm_skew * (Ri)' * Jwi - Ri * rcm_skew * (dRij)' * Jwi - Ri * rcm_skew * (Ri)' * dJwi;
end

function dMj = dMj_Matrix_Internal ( m , dJgi , Jgi , dJwi , Ri , Ig , Jwi , dRij )
    dMj = m * ( (dJgi)' * Jgi + (Jgi)' * dJgi ) + (dJwi)' * Ri * Ig * (Ri)' * Jwi + (Jwi)' * dRij * Ig * (Ri)' * Jwi ...
    + (Jwi)' * Ri * Ig * (dRij)' * Jwi + (Jwi)' * Ri * Ig * (Ri)' * dJwi;
end

function Mp = Mp_Matrix_Internal ( dMj , qp )
    Mp = dMj * qp;
end