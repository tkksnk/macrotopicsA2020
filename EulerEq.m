function  res = EulerEq(cons)
% cを与えたときのオイラー方程式の残差を返す関数

global beta gamma alpha delta capital cfcn0 kgrid

wealth = capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(kgrid(1),kprime);

% 次期の政策関数を線形補間
%cnext = interp1(m.kgrid,cfcn,kprime,'linear','extrap');
% 次期の価値関数をスプライン補間
cnext = interp1(kgrid,cfcn0,kprime,'spline');

%% オイラー方程式
res = (1/cons) - beta*(1/cnext)*(alpha*kprime.^(alpha-1) + (1.-delta));
 
return