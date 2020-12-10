function [coeffs, fun, F] = nep2
%NEP2  3-by-3 basic NEP example.
%  [COEFFS,FUN, F] = nlevp('nep2') generates a 3-by-3
%  nonlinear matrix function, which is equivalent to
%  diag([cos(lambda), sin(lambda), exp(lam) - 7].
%  This problem has the eigenvalues -pi*log(7) and k*pi/2 for any integer k.
%  The matrices are returned in a cell array COEFFS = {A0, A1, A2, A3,
%  A4, A5, A6, A7, A8, A9} and in a function_handle F. The order
%  of the functions of the coefficients is 
%  FUN= {lam, exp(lam), lam*exp(lam), lam*exp(lam)*cos(lam), cos(lam), 
%  lam*cos(lam), sin(lam), lam^2*sin(lam), 1}.
%  This problem has the property nep.

%  Reference: 
%  J. Asakura, T. Sakurai, H. Tadano, T. Ikegami, K. Kimura.
%  A linearization method for nonlinear eigenvalue problems 
%  using a contour integral. Power point prresentation given by Asakura at
%  the 7th International Workshop on Accurate Solution of Eigenvalue Problems
%  Dubrovnik, June 9-12, 2008.
%  http://lavica.fesb.unist.hr/~slap/iwasep7/Talks/iwasep7_asakura.ppt

A0 = [0 0 0; -7 0 -7; 0 0 0];         % z
A1 = [2 0 2; 3 0 3; 1 0 1];           % exp(z)
A2 = [0 0 0; 1 0 1; 0 0 0];           % z*exp(z)
A3 = [0 2 0; 0 3 0; 0 1 0];           % exp(z)*cos(z)
A4 = [0 0 0; 0 1 0; 0 0 0];           % z*exp(z)*cos(z)
A5 = [1 14 0; 0 -21 0; 0 -7 0];       % cos(z)
A6 = [0 0 0; 0 -7 0; 0 0 0];          % z*cos(z)
A7 = [0 -1 0; 0 1 0; 0 0 0];          % sin(z)
A8 = [0 1 0; 0 0 0; 0 0 0];           % z^2*sin(z)
A9 = [-14 0 -14; -21 0 -21; -7 0 -7]; % 1

coeffs = {A0, A1, A2, A3, A4, A5, A6, A7, A8 A9};
fun = @nep2_fun;
F = @(z) [2*exp(z)+cos(z)-14 (z^2-1)*sin(z)+ ...
		2*(exp(z)+7)*cos(z) 2*(exp(z)-7); ...
		(z+3)*(exp(z)-7) sin(z)+ ...
        (z+3)*(exp(z)-7)*cos(z) (z+3)* ...
		(exp(z)-7);exp(z)-7 (exp(z)-7)* ...
        cos(z) exp(z)-7];  
end

function F = nep2_fun(z)
z = z(:);
f0 = z;
f1 = exp(z);
f2 = z.*exp(z);
f3 = exp(z).*cos(z);
f4 = z.*exp(z).*cos(z);
f5 = cos(z);
f6 = z.*cos(z);
f7 = sin(z);
f8 = z.^2.*sin(z);
f9 = ones(size(z));
  
F = [f0 f1 f2 f3 f4 f5 f6 f7 f8 f9];
end
