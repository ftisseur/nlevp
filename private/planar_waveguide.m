function [coeffs,fun,F,xcoeffs] = planar_waveguide
%PLANAR_WAVEGUIDE   Quartic PEP from planar waveguide.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('planar_waveguide') returns the
%  coefficients of a 129-by-129 quartic matrix polynomial
%  lambda^4*A4 + lambda^3*A3 + lambda^2*A2 + lambda*A1 + A0
%  derived from a finite element solution of the propagation constants in a
%  six-layer planar waveguide.  The formulation involves a change of
%  variable, so the actual propagation constants are found by converting
%  the computed eigenvalues; see the reference for details.
%  The matrices are returned in a cell array: COEFFS = {A0,A1,A2,A3,A4}.
%  FUN is a function handle to evaluate the monomials 1,lambda,..,lambda^4
%  and their derivatives.
%  F is the function handle  lambda^4*A4 + lambda^3*A3 + lambda^2*A2 +
%  lambda*A1 + A0.
%  XCOEFFS is the cell {1, L2, 1, L4, 1; A0, R2', A2, R4', A4};
%  with A1 = L1*R1' and A3 = L3*R3', to exploit the low rank of A1 and A3.
%  This problem has the properties pep, real, symmetric, tridiagonal,
%  banded, low-rank. 



%  Reference:
%  D. Stowell and J. Tausch, Variational Formulation for Guided and Leaky
%  Modes in MuliLayer Dielectric Waveguides, Commun. Comput. Phys. 7,2010,
%  pp. 564-579

%  Based on function supplied by David Stowell.

% Frequency = 2*pi/lambda_0 = 2*pi/.6328
w = 9.92918;

% Define the mesh size.
a = 0.0;
b = 2.0;
h = 2^(-6);
n = (b-a)/h;
nl = 6; % Number of layers in the waveguide, including the semi-infinite layers
nlocal = (n/4) -1; % nlocal is the number of interior meshpoints in each layer.

% Call setInfo to set waveguide parameters.
[K,K2,k02,width,deltasq,q] = setInfo(w,nl);

% For Octave compatibility, now assigned at end instead.
% coeffs = {builda0, builda1, builda2, builda3, builda4};


%==========================================================================
% Nested functions.
% These have been un-nested for Octave compatibility.

%=========================================================================
%    function A0 = builda0

       v = 4*ones(n+1,1);
       v(1) = 2;
       v(n+1) = 2;

       u = ones(n,1);
       A0 = diag(v,0) + diag(u,-1) +diag(u,1);
       A0 = (h/6)*(deltasq*deltasq/16)*A0;

%    end
%=========================================================================
%    function A1 = builda1

        u = zeros(n+1,1);
        u(1) = -deltasq/4;
        u(n+1) = deltasq/4;

        A1 = diag(u);

%    end
%=========================================================================
%    function A2 = builda2

       % A2 is built in two steps, A2a and A2b.

       % First A2a.
       v = 2*ones(n+1,1);
       v(1) = 1;
       v(n+1) = 1;
       u = -1*ones(n,1);

       A2a = diag(v,0) + diag(u,-1) +diag(u,1);
       A2a =(1/h)*A2a;

       % A2b
       e = ones(nlocal,1);
       md = zeros(n+1,1);
       supd = zeros(n+1,1);
       subd = zeros(n+1,1);

       md(1) = 2*q(2);
       supd(2) = 1*q(2);
       subd(1) = 1*q(2);

       for k = 1:nl-2

           end_ct = k*(nlocal+1);
           start_ct = end_ct -nlocal+1;

           md(start_ct:end_ct) = 4*q(k+1)*e;
           supd(start_ct +1: end_ct+1) = q(k+1)*e;
           subd(start_ct:end_ct) = q(k+1)*e;

           if k < 4 % Interface points.
               md(end_ct+1) = 4*(q(k+1) + q(k+2))/2.0;
               supd(end_ct+2) = 1*q(k+2);
               subd(end_ct+1) = 1*q(k+2);
           end

       end

       md(n+1) = 2*q(nl-1);
       supd(n+1) = 1*q(nl-1);
       subd(n+1) = 1*q(nl-1);

       A2b = diag(md,0) + diag(subd(1:n),-1) + diag(supd(2:n+1),1);
       A2b = (h/6)*A2b;

       % Now A2.
       A2 = A2a - A2b;

%    end
%=========================================================================
%    function A3 = builda3

        A3 = zeros(n+1,n+1);
        A3(1,1) = 1;
        A3(n+1,n+1) = 1;

%    end
%=========================================================================
%    function A4 = builda4

       v = 4*ones(n+1,1);
       v(1) = 2;
       v(n+1) = 2;

       u = ones(n,1);
       A4 = diag(v,0) + diag(u,-1) +diag(u,1);
       A4 =(h/6)*A4;

%    end
%=========================================================================

coeffs = {A0, A1, A2, A3, A4};
fun = @(lam) nlevp_monomials(lam,4);
F = @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*(coeffs{3}+lam*(coeffs{4} +lam*coeffs{5})));
% building xcoeffs
e1 = zeros(n+1,1); e1(1) = 1;
en = zeros(n+1,1); en(end) = 1;
e1z = e1; e1z(1) = -deltasq/4;
enz = en; enz(end) = deltasq/4;
L2 = [e1z enz];
R2 = [e1 en];
L4 = [e1 en];
R4 = L4;
xcoeffs = {1, L2, 1, L4, 1; A0, R2', A2, R4', A4};
end

function [K,K2,k20,width,deltasq,q] = setInfo(w,nl)
   % Set waveguide parameters used in construction of matrices.
   % For a description of these parameters, see reference above.

    nref = zeros(nl,1);
    K = zeros(nl,1);
    K2 = zeros(nl,1);

    % Refractive indices in each layer.
    nref(1) = 1.5;
    nref(2) = 1.66;
    nref(3) = 1.6;
    nref(4) = 1.53;
    nref(5) = 1.66;
    nref(6) = 1.0;

    K2 = w*w*nref.*nref;

    k20 = K2(1);

    K(1) = 0;
    K(2) = K2(2) - K2(1);
    K(3) = K2(3) - K2(1);
    K(4) = K2(4) - K2(1);
    K(5) = K2(5) - K2(1);

    % Width of each layer.
    width = zeros(nl-2,1);
    width(1) = .5;
    width(2) = .5;
    width(3) = .5;
    width(4) = .5;

    deltasq = K2(1) - K2(nl);

    q = K2 - (K2(1) + K2(nl))/2;
end
