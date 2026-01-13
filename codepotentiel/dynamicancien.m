%% Code Workshop 4 -- TSO Célia et MERBOUCHE Mouloud %%

close all
clc
tic
digits(32)
type = [1;1;0;1];
n = 4;
d = [ 0 ; 0.22 ; 0.24 ; 0 ];
alph = [ 0 ; 0 ; 0 ; 0 ] ;
r = [ 0 ; 0 ; 0 ; 0 ] ;
th = [ 0 ; 0 ; 0 ; 0 ] ;
V_g = [ 0 ; 0 ; 9.81 ] ;

m = sym("m",[ n 1 ],'real') ;
xx= sym("xx",[ n 1 ],'real') ;
yy= sym("yy",[ n 1 ],'real') ;
zz= sym("zz",[ n 1 ],'real') ;
xy= sym("xy",[ n 1 ],'real') ;
xz= sym("xz",[ n 1 ],'real') ;
yz= sym("yz",[ n 1 ],'real') ;

for k = 1 : n
Ig(:,:,k) = [ xx(k) xy(k) xz(k) ; xy(k) yy(k) yz(k) ;
xz(k) yz(k) zz(k) ];
end
x= sym("x",[ n 1 ],'real') ;
y= sym("y",[ n 1 ],'real') ;
z= sym("z",[ n 1 ],'real') ;
for k = 1:n
    rcm(:, k) = [x(k); y(k); z(k)];
end

for k = 1:n
    rcm_skew(:,:,k) = skew(rcm(:, k));
end


Fc = sym('Fc', [n 1], 'real');
Fv = sym('Fv', [n 1], 'real'); 
%p.14
load qf.mat
load dqf.mat
load ddqf.mat
load tauf.mat
load('Data-20251124/q.mat')
load('Data-20251124/tau.mat')

J = 4;
qf   = qf(1:J, :);
dqf  = dqf(1:J, :);
ddqf = ddqf(1:J, :);
tauf = tauf(1:J, :);
q_meas = q_meas(1:J, :);
Tau = Tau(1:J, :);

%p.15
nn = (length(t)-1)/380;
for k = 1 : nn
t_sub(1,k) = t(1,(k-1)*380+1);
q_sub(:,k) = qf(:,(k-1)*380+1);
dq_sub(:,k) = dqf(:,(k-1)*380+1);
ddq_sub(:,k) = ddqf(:,(k-1)*380+1);
tau_sub(:,k) = tauf(:,(k-1)*380+1);
end

tau_id = sym(zeros(n,nn));

for mm = 1:2
%for mm = 1:10 si pas réussi a faire l'optimisation
mm
th(1) = q_sub(1,mm);
th(2) = q_sub(2,mm);
r(3) = q_sub(3,mm);
th(4) = q_sub(4,mm);
qp = dq_sub(:,mm);
qpp = ddq_sub(:,mm);

%p.18
for k = 1:n
        T(:,:,k) = [ cos(th(k))               , - sin(th(k))              , 0             , d(k) ;
                    cos(alph(k)) * sin(th(k)) , cos(alph(k)) * cos(th(k)) , -sin(alph(k)) , - r(k) * sin(alph(k)) ;
                    sin(alph(k)) * sin(th(k)) , sin(alph(k)) * cos(th(k)) , cos(alph(k))  , r(k) * cos(alph(k)) ;
                    0                         , 0                         , 0             , 1 ];
end

Ti(:,:,1) = T(:,:,1);
for k = 2 : n
Ti(:,:,k) = Ti(:,:,k-1) * T(:,:,k);
end

zi(:,1:n) = Ti(1:3,3,1:n);
ti(:,1:n) = Ti(1:3,4,1:n);
Ri(:,:,1:n) = Ti(1:3,1:3,1:n);

%p.20-21
for i = 1 : n
for k = 1 : i
li(:,k,i) = ti(:,i) - ti(:,k);
Jvi(:,k,i) = cross(zi(:,k) , li(:,k,i));
Jwi(:,k,i) = zi(:,k);
end
end

for i = 1 : n
Jgi(:,:,i) = Jvi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i);
end

M = m(1) * (Jgi(:,:,1))' * Jgi(:,:,1) + (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1);
for k = 2 : n
M = M + m(k) * (Jgi(:,:,k))' * Jgi(:,:,k) + (Jwi(:,:,k))' * Ri(:,:,k) * Ig(:,:,k) * (Ri(:,:,k))' *Jwi(:,:,k);
end

%M = m(1) * (Jgi(:,:,1))' * Jgi(:,:,1) + (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1);
%for k = 2 : n
%M = M + m(k) * (Jgi(:,:,k))' * Jgi(:,:,k) + (Jwi(:,:,k))' * Ri(:,:,k) * Ig(:,:,k) * (Ri(:,:,k))' * Jwi(:,:,k);
%end


%p.23-24
for k = 1 : n
if type(k) == 1
dT(:,:,k) = [ - sin(th(k)) , - cos(th(k)) , 0 , 0 ;
cos(alph(k)) * cos(th(k)) , - cos(alph(k)) * sin(th(k)) , 0 , 0 ;
sin(alph(k)) * cos(th(k)) , - sin(alph(k)) * sin(th(k)), 0 , 0 ;
0 , 0 , 0 , 1 ];
elseif type(k) == 0
dT(:,:,k) = [ 0 , 0 , 0 , 0 ;
0 , 0 , 0 , - sin(alph(k)) ;
0 , 0 , 0 , cos(alph(k)) ;
0 , 0 , 0 , 1 ];
end
end


%P.25
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
G(j,1) = m(1) * V_g' * (dtij(:,1,j) + dRij(:,:,1,j) * rcm(:,1));
for i = 2 : n
G(j,1) = G(j,1) + m(i) * V_g' * (dtij(:,i,j) + dRij(:,:,i,j) * rcm(:,i));
end
end

%for j = 1 : n
%G(j,1) = m(1) * V_g' * (dtij(:,1,j) + dRij(:,:,1,j) * rcm(:,1));
%for i = 2 : n
%G(j,1) = G(j,1) + m(i) * V_g' * (dtij(:,i,j) + dRij(:,:,i,j) * rcm(:,i));
%end
%end

for i = 1 : n
for k = 1 : i
for j = 1 : n
dlij(:,k,i,j) = dtij(:,i,j) - dtij(:,k,j);
dJvi(:,k,i,j) = cross(dzij(:,k,j),li(:,k,i)) + cross(zi(:,k),dlij(:,k,i,j));
dJwi(:,k,i,j)  = dzij(:,k,j); 
end
end
end


for i = 1 : n
for j = 1 : n
dJgi(:,:,i,j) = dJvi(:,:,i,j) - dRij(:,:,i,j) * rcm_skew(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (dRij(:,:,i,j))' * Jwi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (Ri(:,:,i))' * dJwi(:,:,i,j);
%dJgi(:,:,i,j) = dJvi(:,:,i,j) - dRij(:,:,i,j) * rcm_skew(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (dRij(:,:,i,j))' * Jwi(:,:,i) - Ri(:,:,i) * rcm_skew(:,:,i) * (Ri(:,:,i))' * dJwi(:,:,i,j);
end
end


for j = 1 : n
    % Initialisation avec le premier corps (i=1) - Page 31
    dMj(:,:,j) = m(1) * ( (dJgi(:,:,1,j))' * Jgi(:,:,1) + (Jgi(:,:,1))' * dJgi(:,:,1,j) ) + ...
                 (dJwi(:,:,1,j))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * Jwi(:,:,1) + ...
                 (Jwi(:,:,1))' * dRij(:,:,1,j) * Ig(:,:,1) * (Ri(:,:,1))'* Jwi(:,:,1) + ...
                 (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (dRij(:,:,1,j))' * Jwi(:,:,1) + ...
                 (Jwi(:,:,1))' * Ri(:,:,1) * Ig(:,:,1) * (Ri(:,:,1))' * dJwi(:,:,1,j);
    
    % Sommation pour les corps suivants (i=2 à n) - Page 32
    for i = 2 : n
        dMj(:,:,j) = (dMj(:,:,j) + m(i) * ( (dJgi(:,:,i,j))' * Jgi(:,:,i) + (Jgi(:,:,i))' * dJgi(:,:,i,j) ) + ...
                     (dJwi(:,:,i,j))' * Ri(:,:,i) * Ig(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) + ...
                     (Jwi(:,:,i))' * dRij(:,:,i,j) * Ig(:,:,i) * (Ri(:,:,i))' * Jwi(:,:,i) + ...
                     (Jwi(:,:,i))' * Ri(:,:,i) * Ig(:,:,i) * (dRij(:,:,i,j))' * Jwi(:,:,i) + ...
                     (Jwi(:,:,i))' * Ri(:,:,i) * Ig(:,:,i) * (Ri(:,:,i))' * dJwi(:,:,i,j));
    end
end 

Mp = simplify(vpa(dMj(:,:,1) * qp(1)));
%Mp = dMj(:,:,1) * qp(1);
for j = 2 : n
    Mp = vpa(Mp + dMj(:,:,j) * qp(j));
    %Mp = Mp + dMj(:,:,j) * qp(j);
end


for j = 1 : n
N2(j,:) = vpa(qp' * dMj(:,:,j));
%N2(j,:) = qp' * dMj(:,:,j);
end
N = vpa(Mp - 1/2 * N2);
%N = Mp - 1/2 * N2;


%tau_id(:,mm) = simplify((M * qpp + N * qp + G - Fv .* qp - Fc .* sign(qp)));
tau_id(:,mm) = simplify(vpa(M * qpp + N * qp + G - Fv .* qp - Fc .* sign(qp)));
end

save tau_id tau_id
save tau_sub.mat tau_sub

toc

%%%%% functions %%%%%
function S=skew(v)
S = [  0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)    0  ];
end 