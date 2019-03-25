function [coeffs,fun,F,xcoeffs] = pillbox_cavity
%PILLBOX_CAVITY  170562-by-170562 NEP from a RF pillbox cavity.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('pillbox_cavity') returns a sparse
%  170562-by-170562 matrices K, M, W1, W2, W3, W4 and W5 for the
%  nonlinear eigenvalue problem defined by T(lam)*x = [K - lam*M +
%  i*(lam-k1)^(1/2)*W1 + i*lam/(lam-k2)^(1/2)*W2 + i*(lam-k3)^(1/2)*W3
%  + i*(lam-k4)^(1/2)*W4 + i*lam/(lam-k4)^(1/2)*W5]*x = 0.
%  COEFFS returns the matrices {K,M,W1,W2,W3,W4,W5} in a cell array.
%  FUN is a function handle to evaluate the functions 1, -lam,
%  i*(lam-k1)^(1/2), i*lam/(lam-k2)^(1/2), i*(lam-k3)^(1/2),
%  i*(lam-k4)^(1/2) and i*lam/(lam-k4)^(1/2), and their derivatives.
%  F is the function handle T(lam).
%  XCOEFFS is the cell {1, 1, L_1, L_2, L_3, L_4, L_5; coeffs{1}, coeffs{2},
%  R_3', R_4', R_5', R_6', R_7'}  that exploits the low rank of Wi =
%  L_1*R_i'.
%  The problem has the property nep, real, symmetryc, sparse, low-rank,
%  banded (31603 bands).

%  Reference: Equation 16 from: 
%  R. Van Beeumen, O. Marques, E.G. Ng, C. Yang, Z. Bai, L. Ge, O. Kononenko,
%  Z. Li, C.-K. Ng, and L. Xiao. Computing resonant modes of accelerator 
%  cavities by solving nonlinear eigenvalue problems via rational
%  approximation, J. Comput. Phys., 374:1031-1043, 2018.

load pillbox_cavity

k1 = 3.625201105315107e+02;
k2 = 4.779508476756487e+02;
k3 = 7.709410421638146e+02;
k4 = 1.581040692232823e+03;

W1 = L3*R3';
W2 = L4*R4';
W3 = L5*R5';
W4 = L6*R6';
W5 = L7*R7';

coeffs = {K, M, W1, W2, W3, W4, W5};
fun = @(lam) pillbox_fun(lam);
F = @(lam) K -lam*M + 1i*(lam -k1)^0.5*W1 +1i*lam/(lam-k2)^0.5*W2...
+ 1i*(lam-k3)^0.5*W3 + 1i*(lam-k4)^0.5*W4 + 1i*lam/(lam-k4)^0.5*W5;

xcoeffs = {1 1 L3 L4 L5 L6 L7; K M R3' R4' R5' R6' R7'};

end

function varargout = pillbox_fun(lam)

lam = lam(:);
n = length(lam);

% cutoff frequencies
k1 = 3.625201105315107e+02;
k2 = 4.779508476756487e+02;
k3 = 7.709410421638146e+02;
k4 = 1.581040692232823e+03;

f1 = 1i*sqrt(lam - k1);
f2 = 1i*lam./sqrt(lam - k2);
f3 = 1i*sqrt(lam - k3);
f4 = 1i*sqrt(lam - k4);
f5 = 1i*lam./sqrt(lam - k4); % this is not a typo.

varargout{1} = [ones(n,1),-lam,f1,f2,f3,f4,f5];
if nargout >= 2

    f1 = 0.5*f4./(lam-k1);
    f21 = 1i*( 0.5./(lam-k2)^(.5));
    f22 = 1i*0.5*k2./(lam-k2)^(1.5);
    f2 = f21 - f22;           % derivative of (lam)/(lam-k2)^.5
    f3 = 0.5*f3./(lam-k3);
    f4 = 0.5*f4./(lam-k4);
    f51 = 1i*( 0.5./(lam-k4)^(.5));
    f52 = 1i*0.5*k2./(lam-k4)^(1.5);
    f5 = f51 - f52;
    varargout{2} = [zeros(n,1),-ones(n,1),f1,f2,f3,f4,f5];
end
for i=2:nargout-1
    f1 = (3/2-i)*f1./(lam-k1);
    f21 = (3/2-i)*f21./(lam-k2);
    f22 = (1/2-i)*f21./(lam-k2);
    f2 = f21 - f22;
    f3 = (3/2-i)*f3./(lam-k3);
    f4 = (3/2-i)*f4./(lam-k4);
    f51 = (3/2-i)*f51./(lam-k4);
    f52 = (1/2-i)*f51./(lam-k4);
    f5 = f51 - f52;
    varargout{i+1} = [zeros(n,2),f1,f2,f3,f4,f5];
end
end
