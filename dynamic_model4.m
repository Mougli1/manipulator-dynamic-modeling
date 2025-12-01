clear
close all
clc

%% === 1. INITIALISATION ===
n = 4; % nombre de ddl

% Paramètres du robot
d = [ 0 ; 0.22 ; 0.24 ; 0 ];
alph = [ 0 ; 0 ; 0 ; 0 ];
r = [ 0 ; 0 ; 0 ; 0 ];
th = [ 0 ; 0 ; 0 ; 0 ];
Vg = [ 0 ; 0 ; 9.81 ];
type = [ 1 ; 1 ; 0 ; 1 ]; % 1=Rotoïde, 0=Prismatique

digits(4);  % Précision pour "vpa"

% Variables symboliques
m = sym("m",[ n 1 ],'real');
xx = sym("xx",[ n 1 ],'real');
yy = sym("yy",[ n 1 ],'real');
zz = sym("zz",[ n 1 ],'real');
xy = sym("xy",[ n 1 ],'real');
xz = sym("xz",[ n 1 ],'real');
yz = sym("yz",[ n 1 ],'real');

% Tenseurs d'inertie (Ig)
for k = 1 : n
    Ig(:,:,k) = [ xx(k) xy(k) xz(k) ; xy(k) yy(k) yz(k) ; xz(k) yz(k) zz(k) ];
end

% Position des centres de masse (rcm)
x = sym("x",[ n 1 ],'real');
y = sym("y",[ n 1 ],'real');
z = sym("z",[ n 1 ],'real');
rcm = [ x' ; y' ; z' ];

% Matrices asymétriques de rcm (nécessite skew.m dans le dossier)
for k = 1 : n
    rcm_skew(:,:,k) = skew(rcm(:,k));
end

Fc = sym("Fc",[ n 1 ],'real');
Fv = sym("Fv",[ n 1 ],'real');

%% === 2. CHARGEMENT DES DONNÉES ===
% Assurez-vous que ces fichiers existent (générés par votre premier script)
% Si vos fichiers sont dans un sous-dossier (ex: Data-2025...), ajoutez le chemin.
if exist('Data-20251124/q.mat', 'file')
    load('Data-20251124/q.mat'); % Charge t et q_meas
else
    load('q.mat'); 
end

load('qf.mat'); 
load('dqf.mat'); 
load('ddqf.mat');
load('tauf.mat');

% Correction de l'erreur "Size inputs must be integers"
% On utilise floor() pour garantir que nn est un entier
nn = floor((length(t)-1)/380); 

% Sous-échantillonnage des données
for k = 1 : nn
    idx = (k-1)*350 + 1; % Indice de sous-échantillonnage
    
    % Vérification pour éviter de dépasser la taille des vecteurs
    if idx > length(t)
        nn = k - 1; % Ajuste nn si on dépasse
        break;
    end
    
    tt(1,k) = t(1, idx);
    % On utilise les données filtrées (qf, dqf...)
    q_sub(:,k) = qf(:, idx);       
    dq_sub(:,k) = dqf(:, idx);
    ddq_sub(:,k) = ddqf(:, idx);
    tau_sub_sampled(:,k) = tauf(:, idx);
end

Tau_id = sym(zeros(n,nn)); % Initialisation de la matrice des couples identifiés

tic
fprintf('Calcul du Modèle Dynamique en cours (%d itérations)...\n', nn);

%% === 3. BOUCLE PRINCIPALE (Dynamic Model) ===
for mm = 1 : 4
    % Affichage de la progression toutes les 10 itérations
    if mod(mm, 10) == 0
        fprintf('Itération %d / %d\n', mm, nn);
    end
    
    % --- Mise à jour des variables d'état (adapté de Loop.m) ---
    th(1) = q_sub(1,mm);
    th(2) = q_sub(2,mm);
    r(3)  = q_sub(3,mm);
    th(4) = q_sub(4,mm);
    
    qp  = dq_sub(:,mm);
    qpp = ddq_sub(:,mm);

    % --- Calcul de T_Matrix ---
    for k = 1 : n
        ct = cos(th(k)); st = sin(th(k));
        ca = cos(alph(k)); sa = sin(alph(k));
        
        T(:,:,k) = [ ct             , -st           , 0          , d(k)             ; 
                     ca * st        , ca * ct       , -sa        , -r(k) * sa       ;
                     sa * st        , sa * ct       , ca         , r(k) * ca        ;             
                     0              , 0             , 0          , 1                ];
    end

    % --- Cinématique Directe (Ti) ---
    Ti(:,:,1) = T(:,:,1);
    for k = 2 : n
        Ti(:,:,k) = Ti(:,:,k-1) * T(:,:,k);
    end

    zi(:,1:n) = Ti(1:3,3,1:n);
    ti(:,1:n) = Ti(1:3,4,1:n);
    Ri(:,:,1:n) = Ti(1:3,1:3,1:n);
 
    % --- Jacobienne (Jacobian_Matrix logic) ---
    for i = 1 : n
        for k = 1 : i
            li(:,k,i) = ti(:,i) - ti(:,k);
            Jvi(:,k,i) = cross(zi(:,k), li(:,k,i));
            Jwi(:,k,i) = zi(:,k);
        end
    end

    % --- Jgi (Jgi_Matrix logic) ---
    for i = 1 : n
        Jgi(:,:,i) = Jvi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i);
    end
    Jgi = simplify(Jgi);

    % --- Matrice d'Inertie M (Inertia_Matrix logic) ---
    % M = m * Jgi' * Jgi + Jwi' * Ri * Ig * Ri' * Jwi
    M = m(1) * (Jgi(:,:,1))' * Jgi(:,:,1) + (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1);
    for k = 2 : n
        M_k = m(k) * (Jgi(:,:,k))' * Jgi(:,:,k) + (Jwi(:,:,k))' * Ri(:,:,k) * Ig(:,:,k) * (Ri(:,:,k))' * Jwi(:,:,k);
        M = M + M_k;
    end
    M = simplify(M);

    % --- dT_Matrix ---
    for k = 1 : n
        if type(k) == 1
            dT(:,:,k) = [ -sin(th(k))             , -cos(th(k))             , 0 , 0 ; 
                          cos(alph(k)) * cos(th(k)) , -cos(alph(k)) * sin(th(k)) , 0 , 0 ;
                          sin(alph(k)) * cos(th(k)) , -sin(alph(k)) * sin(th(k)) , 0 , 0 ;             
                          0                       , 0                       , 0 , 1 ];
        elseif type(k) == 0
            dT(:,:,k) = [ 0 , 0 , 0 , 0 ; 
                          0 , 0 , 0 , -sin(alph(k)) ;
                          0 , 0 , 0 , cos(alph(k)) ;             
                          0 , 0 , 0 , 1 ];
        end
    end

    % --- Calculs dTij, dRij, dzij, dtij ---
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

    % --- Vecteur Gravité G (Gravity logic) ---
    G = sym(zeros(n,1)); 
    for j = 1 : n
        % G = m * Vg' * (dtij + dRij * rcm)
        G(j,1) = m(1) * Vg' * (dtij(:,1,j) + dRij(:,:,1,j) * rcm(:,1));
        for i = 2 : n
            G_term = m(i) * Vg' * (dtij(:,i,j) + dRij(:,:,i,j) * rcm(:,i));
            G(j,1) = G(j,1) + G_term;
        end
    end

    % --- Dérivée de la Jacobienne (dJ_Matrix logic) ---
    for i = 1 : n
        for k = 1 : i
            for j = 1 : n
                dlij(:,k,i,j) = dtij(:,i,j) - dtij(:,k,j);
                % dJvi = cross(dzij,li) + cross(zi,dlij); dJwi = dzij;
                dJvi(:,:,i,j) = cross(dzij(:,k,j), li(:,k,i)) + cross(zi(:,k), dlij(:,k,i,j));
                dJwi(:,:,i,j) = dzij(:,k,j);
            end
        end
    end

    % --- dJgi (dJgi_Matrix logic) ---
    for i = 1 : n
        for j = 1 : n
            dJgi(:,:,i,j) = dJvi(:,:,i,j) ...
                            - dRij(:,:,i,j) * rcm_skew(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) ...
                            - Ri(:,:,i) * rcm_skew(:,:,i) * (dRij(:,:,i,j))' * Jwi(:,:,i) ...
                            - Ri(:,:,i) * rcm_skew(:,:,i) * (Ri(:,:,i))' * dJwi(:,:,i,j);
        end
    end

    % --- dMj (dMj_Matrix logic) ---
    for j = 1 : n
        % Calcul pour k=1
        dMj(:,:,j) = m(1) * ( (dJgi(:,:,1,j))' * Jgi(:,:,1) + (Jgi(:,:,1))' * dJgi(:,:,1,j) ) ...
                     + (dJwi(:,:,1,j))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1) ...
                     + (Jwi(:,:,1))' * dRij(:,:,1,j) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1) ...
                     + (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (dRij(:,:,1,j))' * Jwi(:,:,1) ...
                     + (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * dJwi(:,:,1,j);
        % Calcul pour k=2 à n
        for i = 2 : n
            dMj_term = m(i) * ( (dJgi(:,:,i,j))' * Jgi(:,:,i) + (Jgi(:,:,i))' * dJgi(:,:,i,j) ) ...
                       + (dJwi(:,:,i,j))' * Ri(:,:,i) * Ig(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) ...
                       + (Jwi(:,:,i))' * dRij(:,:,i,j) * Ig(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) ...
                       + (Jwi(:,:,i))' * Ri(:,:,i) * Ig(:,:,i) * (dRij(:,:,i,j))' * Jwi(:,:,i) ...
                       + (Jwi(:,:,i))' * Ri(:,:,i) * Ig(:,:,i) * (Ri(:,:,i))' * dJwi(:,:,i,j);
            dMj(:,:,j) = dMj(:,:,j) + dMj_term;
        end
    end

    % --- Coriolis Mp (Mp_Matrix logic) ---
    Mp = dMj(:,:,1) * qp(1);
    for j = 2 : n
        Mp = Mp + dMj(:,:,j) * qp(j);
    end

    % --- Coriolis N2 et Matrice N ---
    for j = 1 : n
        N2(j,:) = qp' * dMj(:,:,j);
    end
  
    N = Mp - 1/2 * N2;

    % --- Calcul du Couple final Tau_id ---
    % Tau = M*qpp + N*qp + G - Frottements
    Tau_val = M * qpp + N * qp + G - Fv .* qp - Fc .* sign(qp);
    Tau_id(:,mm) = simplify(vpa(Tau_val));

end

%% === 4. SAUVEGARDE ===
save('Tau_id.mat', 'Tau_id');
save('tau_sub_sampled.mat', 'tau_sub_sampled');
toc