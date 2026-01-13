clear, 
close all 
clc 
n = 4;
m = sym("m",[ n 1 ],'real');
x = sym("x",[ n 1 ],'real');
y = sym("y",[ n 1 ],'real');
z = sym("z",[ n 1 ],'real');
xx= sym("xx",[ n 1 ],'real') ;
yy= sym("yy",[ n 1 ],'real') ;
zz= sym("zz",[ n 1 ],'real') ;
xy= sym("xy",[ n 1 ],'real') ;
xz= sym("xz",[ n 1 ],'real') ;
yz= sym("yz",[ n 1 ],'real') ;
Fc = sym("Fc",[ n 1 ],'real');
Fv = sym("Fv",[ n 1 ],'real');

load tau_id.mat;
load tau_sub.mat; 

[ row , col ] = size(tau_sub);
tau_id_conca = sym(zeros(row*col,1));
for k = 1 : col
tau_sub_conca(row*(k-1)+1:row*(k-1)+row,1) =tau_sub(:,k);
tau_id_conca(row*(k-1)+1:row*(k-1)+row,1) = tau_id(:,k);
end

mx = sym("mx",[ n 1 ],'real');
my = sym("my",[ n 1 ],'real');
mz = sym("mz",[ n 1 ],'real');
mxy = sym("mxy",[ n 1 ],'real');
myz = sym("myz",[ n 1 ],'real');
mxz = sym("mxz",[ n 1 ],'real');
mx_sq = sym("mx_sq",[ n 1 ],'real');
my_sq = sym("my_sq",[ n 1 ],'real');
mz_sq = sym("mz_sq",[ n 1 ],'real');


old = [ m .* x .* y ; m .* y .* z ; m .* x .* z ; m .* x.* x ; m .* y .* y ; m .* z .* z ; m .* x ; m .* y ; m.* z ];
new = [ mxy ; myz ; mxz ; mx_sq ; my_sq ; mz_sq ; mx ; my ; mz ];
tau_id_conca_new = subs(tau_id_conca,old,new);

eqns = tau_id_conca_new == tau_sub_conca;

for k = 1 : n
p(18*(k-1)+1:18*(k-1)+18,1) = [ m(k) ; xx(k) ; yy(k) ;
zz(k) ; xy(k) ; xz(k) ; yz(k) ; Fc(k) ; Fv(k) ; mx(k) ;
my(k) ; mz(k) ; mx_sq(k) ; my_sq(k) ; mz_sq(k) ; mxy(k) ; mxz(k) ; myz(k) ];
end

[ D , Y ] = equationsToMatrix ( eqns , p );

D_real = round ( double(D) , 5 );
Y_real = double ( Y );

%%

% Find the indices of the columns in matrix D that are
%filled with zeros :
cols_with_all_zeros = find ( all ( D_real == 0 ) );
% Eliminate the parameters that do not affect the
%robot’s dynamics, i.e., those that correspond to zero
%columns in the observation matrix :
p ( cols_with_all_zeros ) = [ ];
% Remove column if the entire column is zero :
D_real ( : , all ( D_real == 0 ) ) = [ ];



% Matrix rank :
r = rank(D_real);
% QR Decomposition :
[ Q, R, E ] = qr(D_real, 'vector');
% Find the independent columns :
idx_ind = E(1:r);
% Reduced observation matrix ’Db’ :
Db = D_real(:, idx_ind);
% Basic parameter vector ’pb’ :
Kd = R(1:r,1:r) \ R(1:r,r+1:end);
pb = p(idx_ind) + Kd * p(E(r+1:end));


% Rounding of the basic parameter vector :
tol = 1e-4;
for k = 1 : length(pb)
%Separate the coefficients and parameters
[coeffsExpr, unknown] = coeffs(pb(k));
%Cancel the coefficients that are less than "tol" :
coeffsExpr(abs(coeffsExpr) < tol) = 0;
%New rounded parameter vector :
pb(k) = simplify(sum(coeffsExpr .* unknown));
end
%%

fun = @(x)(sqrt((Y_real - Db * x)' * (Y_real - Db * x)));
x0 = zeros(length(pb),1);
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e5);
sol = fmincon(fun, x0, [], [], [], [], [], [], [], options);
save('sol.mat','sol');
save('id_struct.mat','idx_ind','p','cols_with_all_zeros');
