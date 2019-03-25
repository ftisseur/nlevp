function [coeffs,fun,F,xcoeffs,sol] = fiber
%FIBER  NEP from fiber optic design.
%  [COEFFS,FUN,F,XCOEFFS,SOL] = nlevp('fiber') returns the 2400-by-2400
%  matrices A, I and B in the nonlinear eigenvalue problem
%  T(lambda)*x = (A - lambda*I + s(lambda)*B)x = 0,
%  where A is symmetric tridiagonal and B is zero except for
%  (2400,2400) element equal to 1.
%  The matrices are returned in a cell array: COEFFS = {A, I, B}.
%  FUN is a function handle to evaluate the functions 1, -lambda, s(lambda)
%  and their first two derivatives.
%  F is the function handle T(lambda) = A - lambda*I +s(lambda)*B. 
%  XCOEFFS is the cell {1, 1, en; A, I, en'}, with B = en*en' to exploit
%  the low rank of B.
%  Computed approximations to the smallest positive eigenvalue and corresponding
%  eigenvector are returned as SOL.EVAL and SOL.EVEC.
%
%  This problem has the properties nep, sparse, solution, tridiagonal,
%  banded, low-rank. 

%  References:
%  X. Huang, Z. Bai, and Y. Su.
%    Nonlinear rank-one modification of the symmetric eigenvalue problem.
%    J. Comput. Math.}, 28(2):218-234, 2010.
%  L. Kaufman. Eigenvalue problems in fiber optic design.
%    SIAM J. Matrix Anal. Appl., 28(1):105-117, 2006.

alpha = 25; gamma = 0.003; ell = 1.1; delta = 0.01;
eta_cl = 1.4969; k_cl = 2*pi*eta_cl/ell;
n_c = 400; n = 6*n_c;
r = (1:n+1)*delta; m = 1;

inc = (1:n_c)';
in = (n_c+1:n-1)';

% The following tree variables C_r, eta_r, and k_r were nested functions in
% previous versions; converted in order to supported Octave.
C_r = sqrt( (1 - 2*gamma*(inc/n_c).^alpha) / (1 - 2*gamma) ) - 1;
eta_r = eta_cl + 1.4201*C_r;
k_r = 2*pi*eta_r/ell;

e = ones(n_c,1);
%y1 = -2*e - m^2*(e ./ inc.^2) + delta^2*(k(r(1:n_c)).^2 - k_cl^2);
y1 = -2*e - m^2*(e ./ inc.^2) + delta^2*(k_r(1:n_c).^2 - k_cl^2);
e = ones(size(in));
y2 = -2*e - m^2*(e ./ in.^2);
y = [y1; y2; -1 + 1/(2*n) - m^2/n^2];
i = 1:n-1;
z = (i+0.5) ./ sqrt( i.*(i+1) );
if nlevp_isoctave
   % Octave has no gallery function
   A = spdiags([ [z(:);0],y(:),[0;z(:)] ], -1:1, n, n);
else
   A = gallery('tridiag',z,y,z);
end

B = sparse(n,n); B(n,n) = 1;
I = speye(n);

coeffs = {A,I,B};
fun = @fiber_fun2;
s = @(lam) bkr1(1, n^2*lam)*(n+0.5)/n^2;
F = @(lam) A -lam.*I + s(lam)*B;

en = sparse(zeros(n,1)); en(n) = 1;

xcoeffs1 = {1,  1, en};
xcoeffs2 = {coeffs{1}, coeffs{2}, en'};
xcoeffs = {xcoeffs1{:}; xcoeffs2{:}};

if nargout >= 5
   load fiber.mat  % Loads sol.
end



%      %%%%%%%%%%%%%%%%%%%
%      % NB: no r dependence in next 3 functions, but this matches
%      %     notation in Huang, Bai and Su.
%      % Nested functions
%      function f = eta(r)
%      f = eta_cl + 1.4201*C(r);
%      end
%
%      function f = C(r)
%      f = sqrt( (1 - 2*gamma*(inc/n_c).^alpha) / (1 - 2*gamma) ) - 1;
%      end
%
%      function f = k(r)
%      f = 2*pi*eta(r)/ell;
%      end

end

function [f,fp,fpp] = fiber_fun2(lam)
lam = lam(:);
n = length(lam);
[s,ds,d2s] = fiber_fun(lam);
f = [ones(n,1),-lam,s];
fp = [zeros(n,1),-ones(n,1),ds];
fpp = [zeros(n,2),d2s];
end

function [s,ds,d2s] = fiber_fun(lambda)
%FIBER_FUN  Scalar function associated with FIBER problem.
%  [s,ds,d2s] = NLEVP('fiber_fun',lambda) evaluate the scalar function
%  s(lambda) arising in the NLEVP('fiber') problem and its first (ds)
%  and second (d2s) derivatives.

%  s(lambda) for lambda >= 0 with mode number m and grid number L is
%
%                  L+0.5                     K'(m,sqrt(lambda)L)
%   s(lambda) =  ------ * sqrt(lambda)L *  ------------------- ,
%                   L^2                      K(m,sqrt(lambda)L)
%
%   where K(m,z) = besselk(m,z) is the m-th order modified Bessel
%   function of the 2nd kind, K' is its derivative to z.

%   Ref: [1] Theory of Bessel Functions, G.N.WATSON, 2nd ed., 1944.

%   Based on code by Xin Huang, Fudan University, China, xinhuang@fudan.edu.cn
%   Revision date: 2010/6/1

n = length(lambda);

s = zeros(size(lambda));
ds = s;
d2s = s;

m = 1;
L = 2400;

for i = 1:n
    [s(i), ds(i), d2s(i)] = bkr(m,lambda(i)*L^2);
    %function bkr(m,x) = sqrt(x)*K'(m,sqrt(x))/K(m,sqrt(x))
end

% --- multiply the factor ---
s = (L+.5)/L^2*s;
ds = (L+.5)*ds;
d2s = (L+.5)*L^2*d2s;

end

% ----------------------------------------------------------------

function [s,ds,d2s] = bkr(m,x)
% bkr(m,x) = sqrt(x)*K'(m,sqrt(x))/K(m,sqrt(x))
% The derivative of the BesselK function is computed using the following
% recurrence property [1,p.79]:
%  K'(m,z) = - K(m-1,z) - m*K(m,z)/z = m*K(m,z)/z - K(m+1,z).
% We only use k0 = K(m-1,z), k1 = K(m,z), with sqrtx = sqrt(x)=sqrt(lambda)*L.

if abs(x)< 1e-30 %lower bound for input x
    x = 1e-30;
end

sqrtx = sqrt(x);
k0 = besselk(m-1,sqrtx);
k1 = besselk(m,sqrtx);

a = k0/k1;

if (k0==0) || (k1==0) %underflow
    %truncated to the 2nd term
    %a = (8*sqrtx+(4*(m-1)^2-1))/(8*sqrtx+(4*m^2-1));

    %truncated to the 3rd term
    a = (128*x+16*sqrtx*(4*(m-1)^2-1)+(4*(m-1)^2-1)*(4*(m-1)^2-9))...
        /(128*x+16*sqrtx*(4*m^2-1)+(4*m^2-1)*(4*m^2-9));

    %trancated to the 5th term
    %m1=m-1;
    %a = (1+(4*m1*m1-1)/8/sqrtx*(1+(4*m1*m1-9)/8/sqrtx*(1/2+(4*m1*m1-25)/8/sqrtx*(1/3+(4*m1*m1-49)/8/sqrtx/4)))) ...
    %    /...
    %    (1+(4*m*m-1)/8/sqrtx*(1+(4*m*m-9)/8/sqrtx*(1/2+(4*m*m-25)/8/sqrtx*(1/3+(4*m*m-49)/8/sqrtx/4))));
end
%When x is large enough, underflow occurs for k0 and k1. From the
%asymptotic expansion [1,p.202]: K(m,z)~(pi/2z)^(1/2)exp(-z)[1+(4m^2-1)/(8z)
%+(4m^2-1)(4m^2-3^2)/(2*(8z)^2)+...], we truncate it to the 3rd term and
%let a = (128z^2+16z*(4*(m-1)^2-1)+(4*(m-1)^2-1)*(4*(m-1)^2-9))/...
%(128*z^2+16z*(4*m^2-1)+(4*m^2-1)*(4*m^2-9)).

b = (a*a+(2*m-1)*a/sqrtx-1)/2/sqrtx;    %b=da/dx

s = -(sqrtx*a+m);
ds = -(a*a+2*m*a/sqrtx-1)/2;
d2s = -a*b+m/2*x^(-3/2)*a-m*b/sqrtx;

end


function s = bkr1(m,x)
% bkr1(m,x) = sqrt(x)*K'(m,sqrt(x))/K(m,sqrt(x)) is a copy of bkr that does
% not compute the derivatives.
% We only use k0 = K(m-1,z), k1 = K(m,z), with sqrtx = sqrt(x)=sqrt(lambda)*L.

if abs(x)< 1e-30 %lower bound for input x
    x = 1e-30;
end

sqrtx = sqrt(x);
k0 = besselk(m-1,sqrtx);
k1 = besselk(m,sqrtx);

a = k0/k1;

if (k0==0) || (k1==0) %underflow
    %truncated to the 2nd term
    %a = (8*sqrtx+(4*(m-1)^2-1))/(8*sqrtx+(4*m^2-1));

    %truncated to the 3rd term
    a = (128*x+16*sqrtx*(4*(m-1)^2-1)+(4*(m-1)^2-1)*(4*(m-1)^2-9))...
        /(128*x+16*sqrtx*(4*m^2-1)+(4*m^2-1)*(4*m^2-9));

    %trancated to the 5th term
    %m1=m-1;
    %a = (1+(4*m1*m1-1)/8/sqrtx*(1+(4*m1*m1-9)/8/sqrtx*(1/2+(4*m1*m1-25)/8/sqrtx*(1/3+(4*m1*m1-49)/8/sqrtx/4)))) ...
    %    /...
    %    (1+(4*m*m-1)/8/sqrtx*(1+(4*m*m-9)/8/sqrtx*(1/2+(4*m*m-25)/8/sqrtx*(1/3+(4*m*m-49)/8/sqrtx/4))));
end
%When x is large enough, underflow occurs for k0 and k1. From the
%asymptotic expansion [1,p.202]: K(m,z)~(pi/2z)^(1/2)exp(-z)[1+(4m^2-1)/(8z)
%+(4m^2-1)(4m^2-3^2)/(2*(8z)^2)+...], we truncate it to the 3rd term and
%let a = (128z^2+16z*(4*(m-1)^2-1)+(4*(m-1)^2-1)*(4*(m-1)^2-9))/...
%(128*z^2+16z*(4*m^2-1)+(4*m^2-1)*(4*m^2-9)).

s = -(sqrtx*a+m);

end

