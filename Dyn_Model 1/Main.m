clear
close all
clc

n = 4; % nombre de ddl

d = [ 0 ; 0.22 ; 0.24 ; 0 ];
alph = [ 0 ; 0 ; 0 ; 0 ];
r = [ 0 ; 0 ; 0 ; 0 ];
th = [ 0 ; 0 ; 0 ; 0 ];
Vg = [ 0 ; 0 ; 9.81 ];
digits(4);  % Précision après la virgule (pour la fonction "vpa")
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

% Data loading :
load q.mat
load qf.mat 
load dqf.mat 
load ddqf.mat
load tauf.mat
nn = (length(t)-1)/380;

% Data sub-sampling :
for k = 1 : nn
    tt(1,k) = t(1,(k-1)*350+1);
    q(:,k) = qf(:,(k-1)*350+1);
    dq(:,k) = dqf(:,(k-1)*350+1);
    ddq(:,k) = ddqf(:,(k-1)*350+1);
    tau_sub_sampled(:,k) = tauf(:,(k-1)*350+1);
end

Tau_id = sym(zeros(n,nn));

tic

% Computing Dynamic Model : 
for mm = 1 : nn
    mm
    Loop
end
save Tau_id.mat Tau_id
save tau_sub_sampled.mat tau_sub_sampled
toc

