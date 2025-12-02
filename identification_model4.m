clear, close all 
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
my(k) ; mz(k) ; mx_sq(k) ; my_sq(k) ; mz_sq(k) ; mxy(k)
; mxz(k) ; myz(k) ];
end

[ D , Y ] = equationsToMatrix ( eqns , p );

D_real = round ( double(D) , 5 );
Y_real = double ( Y );
cols_with_all_zeros = find ( all ( D_real == 0 ) );
p ( cols_with_all_zeros ) = [ ];
D_real ( : , all ( D_real == 0 ) ) = [ ];


[R,idx] = rref(D_real);


s = 0;
for k = 1 : length(idx) - 1
dif = idx(k+1) - idx(k);
if dif ~= 1
vd(s+1:s+dif-1,1) = idx ( k ) + 1 : idx ( k ) + dif - 1;
s = length(vd);
end
end



[ rD , cD ] = size( D_real );
s = 1;
for k = 1 : length(vd)
for mm = 1 : length(idx)
[ roww_1 , coll_1 , D_real_nz_1 ] =find(D_real(:,vd(k)));
[ roww_2 , coll_2 , D_real_nz_2 ] = find(D_real(:,idx(mm)));
if length(roww_1) == length(roww_2)
if roww_1 == roww_2
K = D_real_nz_1 ./ D_real_nz_2;
if ((K(1)*ones(length(roww_1),1)) == K)
DV(s,k) = idx(mm);
s = s + 1;
end
end
end
end
s = 1;
end

for k = 1 : length(vd)
p(DV(k)) = p(DV(k)) + K(1) * p(vd(k));
p(vd(k)) = [ ];
D_real(:,vd(k)-k+1) = [ ];
end

fun = @(x)( sqrt( (Y_real - D_real * x)' * (Y_real - D_real * x) ) );
x0 = zeros ( length(p),1);
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e5);
sol = fmincon(fun, x0, [], [], [], [], [], [], [], options);