%%Standalone mais avec idx etape d'avant

%% validation.m — Workshop 4 DOF - TSO Célia & MERBOUCHE Mouloud (standalone, correct)

clear; close all; clc
digits(32);

% ===== 0) Charger sol + structure d’identification =====
load('sol.mat','sol');
load('id_struct.mat','idx_ind','p');   % idx_ind et p viennent du dataset d’identification

% ===== 1) Charger Validation.mat =====
load('validation.mat');  % q_meas_validation, tau_meas_validation, t

q4 = q_meas_validation(1:4, :);

% --- IMPORTANT: ici tau_meas_validation = COURANTS (selon le prof) ---
I4 = tau_meas_validation(1:4, :);      % courants moteurs

Kc = [0.38; 0.38; 0.22; 0.21];         % Nm/A
Gr = [120; 160; 120; 100];             % gear ratios

% Si jamais I est en mA, décommente:
% I4 = I4 / 1000;

Tau4 = (Kc .* Gr) .* I4;               % couples équivalents articulations

%% 2) Data processing (comme data_processing)
q_conv = q4;
q_conv([1 2 4], :) = deg2rad(q4([1 2 4], :));
q_conv(3, :)       = q4(3, :) / 1000;

fc = 10;
fe = 1000;
Wn = 2 * fc / fe;
[b, a] = butter(1, Wn);

[n_joints, n_samples] = size(q_conv);

qf   = zeros(n_joints, n_samples);
dqf  = zeros(n_joints, n_samples);
ddqf = zeros(n_joints, n_samples);
tauf = zeros(n_joints, n_samples);

for i = 1:n_joints
    qf(i, :)   = filtfilt(b, a, q_conv(i, :));
    tauf(i, :) = filtfilt(b, a, Tau4(i, :));

    dq_raw     = derive(t, qf(i, :));
    dqf(i, :)  = filtfilt(b, a, dq_raw);

    ddq_raw    = derive(t, dqf(i, :));
    ddqf(i, :) = filtfilt(b, a, ddq_raw);
end

%% 3) Modèle dynamique (comme script 2)
n = 4;

type = [1;1;0;1];
d    = [0; 0.22; 0.24; 0];
alph = [0; 0; 0; 0];
r    = [0; 0; 0; 0];
th   = [0; 0; 0; 0];

V_g = [0; 0; 9.81];

m  = sym("m",[n 1],'real');
xx = sym("xx",[n 1],'real'); yy = sym("yy",[n 1],'real'); zz = sym("zz",[n 1],'real');
xy = sym("xy",[n 1],'real'); xz = sym("xz",[n 1],'real'); yz = sym("yz",[n 1],'real');

Ig = sym(zeros(3,3,n));
for k = 1:n
    Ig(:,:,k) = [xx(k) xy(k) xz(k);
                 xy(k) yy(k) yz(k);
                 xz(k) yz(k) zz(k)];
end

x = sym("x",[n 1],'real'); y = sym("y",[n 1],'real'); z = sym("z",[n 1],'real');
rcm = sym(zeros(3,n));
rcm_skew = sym(zeros(3,3,n));
for k = 1:n
    rcm(:,k) = [x(k); y(k); z(k)];
    rcm_skew(:,:,k) = skew(rcm(:,k));
end

Fc = sym('Fc', [n 1], 'real');
Fv = sym('Fv', [n 1], 'real');

%% Sous-échantillonnage (comme script 2)
nn = (length(t)-1)/380;   % supposé entier dans le dataset

t_sub   = zeros(1, nn);
q_sub   = zeros(n, nn);
dq_sub  = zeros(n, nn);
ddq_sub = zeros(n, nn);
tau_sub = zeros(n, nn);

for k = 1:nn
    idx = (k-1)*380 + 1;
    t_sub(1,k)    = t(1,idx);
    q_sub(:,k)    = qf(:,idx);
    dq_sub(:,k)   = dqf(:,idx);
    ddq_sub(:,k)  = ddqf(:,idx);
    tau_sub(:,k)  = tauf(:,idx);
end

tau_id = sym(zeros(n,nn));

for mm = 1:nn
    mm
    th(1) = q_sub(1,mm);
    th(2) = q_sub(2,mm);
    r(3)  = q_sub(3,mm);
    th(4) = q_sub(4,mm);

    qp  = dq_sub(:,mm);
    qpp = ddq_sub(:,mm);

    % --- T, Ti, zi, ti, Ri ---
    T = sym(zeros(4,4,n));
    for k = 1:n
        T(:,:,k) = [ cos(th(k))               , - sin(th(k))              , 0             , d(k) ;
                     cos(alph(k))*sin(th(k))  ,  cos(alph(k))*cos(th(k))  , -sin(alph(k)) , -r(k)*sin(alph(k)) ;
                     sin(alph(k))*sin(th(k))  ,  sin(alph(k))*cos(th(k))  ,  cos(alph(k)) ,  r(k)*cos(alph(k)) ;
                     0                        ,  0                        ,  0            , 1 ];
    end

    Ti = sym(zeros(4,4,n));
    Ti(:,:,1) = T(:,:,1);
    for k = 2:n
        Ti(:,:,k) = Ti(:,:,k-1)*T(:,:,k);
    end

    zi = sym(zeros(3,n));
    ti = sym(zeros(3,n));
    Ri = sym(zeros(3,3,n));
    for k = 1:n
        zi(:,k)   = Ti(1:3,3,k);
        ti(:,k)   = Ti(1:3,4,k);
        Ri(:,:,k) = Ti(1:3,1:3,k);
    end

    % --- Jacobians ---
    Jvi = sym(zeros(3,n,n));
    Jwi = sym(zeros(3,n,n));
    li  = sym(zeros(3,n,n));
    for i = 1:n
        for k = 1:i
            li(:,k,i)  = ti(:,i) - ti(:,k);
            Jvi(:,k,i) = cross(zi(:,k), li(:,k,i));
            Jwi(:,k,i) = zi(:,k);
        end
    end

    Jgi = sym(zeros(3,n,n));
    for i = 1:n
        Jgi(:,:,i) = Jvi(:,:,i) - Ri(:,:,i)*rcm_skew(:,:,i)*(Ri(:,:,i))'*Jwi(:,:,i);
    end

    % --- Mass matrix M ---
    M = m(1)*(Jgi(:,:,1))'*Jgi(:,:,1) + (Jwi(:,:,1))'*Ri(:,:,1)*Ig(:,:,1)*(Ri(:,:,1))'*Jwi(:,:,1);
    for k = 2:n
        M = M + m(k)*(Jgi(:,:,k))'*Jgi(:,:,k) + (Jwi(:,:,k))'*Ri(:,:,k)*Ig(:,:,k)*(Ri(:,:,k))'*Jwi(:,:,k);
    end

    % --- dT ---
    dT = sym(zeros(4,4,n));
    for k = 1:n
        if type(k) == 1
            dT(:,:,k) = [ -sin(th(k)) , -cos(th(k)) , 0 , 0 ;
                           cos(alph(k))*cos(th(k)) , -cos(alph(k))*sin(th(k)) , 0 , 0 ;
                           sin(alph(k))*cos(th(k)) , -sin(alph(k))*sin(th(k)) , 0 , 0 ;
                           0 , 0 , 0 , 1 ];
        else
            dT(:,:,k) = [ 0 , 0 , 0 , 0 ;
                          0 , 0 , 0 , -sin(alph(k)) ;
                          0 , 0 , 0 ,  cos(alph(k)) ;
                          0 , 0 , 0 , 1 ];
        end
    end

    % --- dTij, dRij, dzij, dtij ---
    dTij = sym(zeros(4,4,n,n));
    dRij = sym(zeros(3,3,n,n));
    dzij = sym(zeros(3,n,n));
    dtij = sym(zeros(3,n,n));

    for i = 1:n
        for j = 1:n
            if j > i
                dTij(:,:,i,j) = zeros(4);
            elseif j == i
                if j == 1
                    dTij(:,:,i,j) = dT(:,:,j);
                else
                    dTij(:,:,i,j) = Ti(:,:,j-1)*dT(:,:,j);
                end
            else
                Tij = eye(4);
                for k = j+1:i
                    Tij = Tij*Ti(:,:,k);
                end
                if j == 1
                    dTij(:,:,i,j) = dT(:,:,j)*Tij;
                else
                    dTij(:,:,i,j) = Ti(:,:,j-1)*dT(:,:,j)*Tij;
                end
            end
            dRij(:,:,i,j) = dTij(1:3,1:3,i,j);
            dzij(:,i,j)   = dTij(1:3,3,i,j);
            dtij(:,i,j)   = dTij(1:3,4,i,j);
        end
    end

    % --- Gravity vector G ---
    G = sym(zeros(n,1));
    for j = 1:n
        G(j,1) = m(1)*V_g'*(dtij(:,1,j) + dRij(:,:,1,j)*rcm(:,1));
        for i = 2:n
            G(j,1) = G(j,1) + m(i)*V_g'*(dtij(:,i,j) + dRij(:,:,i,j)*rcm(:,i));
        end
    end

    % --- dJ terms ---
    dlij = sym(zeros(3,n,n,n));
    dJvi = sym(zeros(3,n,n,n));
    dJwi = sym(zeros(3,n,n,n));

    for i = 1:n
        for k = 1:i
            for j = 1:n
                dlij(:,k,i,j) = dtij(:,i,j) - dtij(:,k,j);
                dJvi(:,k,i,j) = cross(dzij(:,k,j), li(:,k,i)) + cross(zi(:,k), dlij(:,k,i,j));
                dJwi(:,k,i,j) = dzij(:,k,j);
            end
        end
    end

    dJgi = sym(zeros(3,n,n,n));
    for i = 1:n
        for j = 1:n
            dJgi(:,:,i,j) = dJvi(:,:,i,j) ...
                          - dRij(:,:,i,j)*rcm_skew(:,:,i)*(Ri(:,:,i))'*Jwi(:,:,i) ...
                          - Ri(:,:,i)*rcm_skew(:,:,i)*(dRij(:,:,i,j))'*Jwi(:,:,i) ...
                          - Ri(:,:,i)*rcm_skew(:,:,i)*(Ri(:,:,i))'*dJwi(:,:,i,j);
        end
    end

    % --- dM ---
    dMj = sym(zeros(n,n,n));
    for j = 1:n
        dMj(:,:,j) = m(1)*((dJgi(:,:,1,j))'*Jgi(:,:,1) + (Jgi(:,:,1))'*dJgi(:,:,1,j)) + ...
                     (dJwi(:,:,1,j))'*Ri(:,:,1)*Ig(:,:,1)*(Ri(:,:,1))'*Jwi(:,:,1) + ...
                     (Jwi(:,:,1))'*dRij(:,:,1,j)*Ig(:,:,1)*(Ri(:,:,1))'*Jwi(:,:,1) + ...
                     (Jwi(:,:,1))'*Ri(:,:,1)*Ig(:,:,1)*(dRij(:,:,1,j))'*Jwi(:,:,1) + ...
                     (Jwi(:,:,1))'*Ri(:,:,1)*Ig(:,:,1)*(Ri(:,:,1))'*dJwi(:,:,1,j);
        for i = 2:n
            dMj(:,:,j) = dMj(:,:,j) + m(i)*((dJgi(:,:,i,j))'*Jgi(:,:,i) + (Jgi(:,:,i))'*dJgi(:,:,i,j)) + ...
                         (dJwi(:,:,i,j))'*Ri(:,:,i)*Ig(:,:,i)*(Ri(:,:,i))'*Jwi(:,:,i) + ...
                         (Jwi(:,:,i))'*dRij(:,:,i,j)*Ig(:,:,i)*(Ri(:,:,i))'*Jwi(:,:,i) + ...
                         (Jwi(:,:,i))'*Ri(:,:,i)*Ig(:,:,i)*(dRij(:,:,i,j))'*Jwi(:,:,i) + ...
                         (Jwi(:,:,i))'*Ri(:,:,i)*Ig(:,:,i)*(Ri(:,:,i))'*dJwi(:,:,i,j);
        end
    end

    Mp = dMj(:,:,1)*qp(1);
    for j = 2:n
        Mp = Mp + dMj(:,:,j)*qp(j);
    end

    N2 = sym(zeros(n,n));
    for j = 1:n
        N2(j,:) = qp' * dMj(:,:,j);
    end
    N = Mp - 1/2*N2;

    tau_id(:,mm) = simplify(vpa(M*qpp + N*qp + G - Fv.*qp - Fc.*sign(qp)));
end

%% 4) Construire Db_val AVEC p et idx_ind de l’identification, puis tau_c
[row, col] = size(tau_sub);

tau_id_conca  = sym(zeros(row*col,1));
tau_sub_conca = zeros(row*col,1);
for k = 1:col
    tau_sub_conca(row*(k-1)+1:row*k,1) = tau_sub(:,k);
    tau_id_conca(row*(k-1)+1:row*k,1)  = tau_id(:,k);
end

% changement de variables (comme script 3)
mx    = sym("mx",[n 1],'real'); my = sym("my",[n 1],'real'); mz = sym("mz",[n 1],'real');
mxy   = sym("mxy",[n 1],'real'); myz = sym("myz",[n 1],'real'); mxz = sym("mxz",[n 1],'real');
mx_sq = sym("mx_sq",[n 1],'real'); my_sq = sym("my_sq",[n 1],'real'); mz_sq = sym("mz_sq",[n 1],'real');

old = [ m.*x.*y ; m.*y.*z ; m.*x.*z ; m.*x.*x ; m.*y.*y ; m.*z.*z ; m.*x ; m.*y ; m.*z ];
new = [ mxy ; myz ; mxz ; mx_sq ; my_sq ; mz_sq ; mx ; my ; mz ];

tau_id_conca_new = subs(tau_id_conca, old, new);
eqns = tau_id_conca_new == tau_sub_conca;

[D_val, ~] = equationsToMatrix(eqns, p);
D_real_val = round(double(D_val), 5);

Db_val = D_real_val(:, idx_ind);

% cohérence dimensions
assert(size(Db_val,2) == length(sol), "Db_val et sol ne sont pas compatibles (idx_ind / sol mismatch).");

tau_c_conca = Db_val * sol;
tau_c = reshape(tau_c_conca, n, nn);

%% 5) Tracés
figure('Name','Validation: mesuré (filtré) vs calculé (Db_val*sol)');
for j = 1:n
    subplot(2,2,j);
    plot(t_sub, tau_sub(j,:), 'LineWidth', 1.2); hold on
    plot(t_sub, tau_c(j,:),   'LineWidth', 1.2);
    grid on
    xlabel('t (s)');
    ylabel('Torque');
    title(['Joint ', num2str(j)]);
    legend('Mesuré (filtré)', 'Calculé', 'Location','best');
end

%%%% functions %%%%%
function df = derive(t, f)
    df = [ f(1,2) - f(1,1), ...
           (f(1,3:end) - f(1,1:end-2)) / 2, ...
           f(1,end) - f(1,end-1) ] / (t(1,2) - t(1,1));
end

function S = skew(v)
    S = [  0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0 ];
end
