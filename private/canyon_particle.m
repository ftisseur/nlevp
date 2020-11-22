function [coeffs, fun, F, xcoeffs] = canyon_particle(stepSize, w1, w2, l, m, depth)
%CANYON_PARTICLE NEP  from the Schroedinger equation on a canyon-like shape.
%  [COEFFS,FUN,F,XCOEFFS] = nlevp('canyon_particle',stepSize,w1,w2,l,m,depth)
%  generates a N-by-N (N is 16281 by default) NEP from the finite element
%  model of the Schrödinger equation on a canyon-like shape. It has the form 
%  F(lam) = H -lam*eye - sum_{k=1}^nz*exp(1i*sqrt(m*(lam-a_k)))*L_k*U_k',
%  where nz is (by default) 81, H is a sparse N-by-N matrix, and L_k and
%  U_k are N-by-2 sparse matrices.
%  The user may set the stepSize (default = 0.05) of the finite element
%  model. The stepSize influences nz, which is given by floor(4/stepSize)+1
%  and N, which is given by nz*nx, where nx = floor(10/stepSize)+1.
%  They may also set the widths of the canyon, w1 < w2 < 2, where the
%  default values are 1 and 1.1, and the length of the canyon l < 10, where
%  the default is 4.
%  Finally, the user may also set the mass of the particle (default 0.2)
%  and the depth of the canyon (default 3).
%  COEFFS is the cell of the coefficients of length nz, COEFFS = {H, eye,
%  L_1*U_1', ..., L_nz*U_nz'}.
%  FUN is a function handle to evaluate the functions 1,-lambda, and the nz
%  functions -exp(1i*sqrt(m*(lam-a_k))).
%  F is the function handle F(lam).
%  XCOEFFS returns the cell {1, 1, L_1, ... L_nz; H, eye, U_1', ..., U_nz'}
%  to exploit the low rank of the problem.
%  The problem has the properties nep, sparse, scalable,
%  parameter-dependent, banded (80 bands), low-rank.
 
%  References:
%  S. G�ttel, R. Van Beeumen, K. Meerbergen, and W. Michiels. 
%  NLEIGS: a class of fully rational Krylov methods for nonlinear eigenvalue 
%  problems. SIAM J. Sci. Comput., 36(6):A2842�A2864, 2014.
%
%  W. Vandenberghe, M. Fischetti, R. Van Beeumen, K. Meerbergen, W. Michiels, 
%  and C. Effenberger. Determining bound states in a semiconductor device with 
%  contacts using a nonlinear eigenvalue solver. 
%  Journal of Computational Electronics, 13(3):753�762, 2014.
%
%  A part of this code is an adaptation from "particle_init.m" by William 
%  Vandenberghe, which can be found in NLEIGS v0.5 (23/01/19).

% Constants
meter = 1/5.2917725e-11;
nm = 1e-9*meter;
eV = 1/13.6;

% parameters
xmax = 5;
zmax = 2;
interval = 2;

% Setting the default inputs if not given
if nargin < 1 || isempty(stepSize)
    stepSize = 0.05;
end
if nargin < 2 || isempty(w1)
   w1 = 1;
end
if nargin < 3 || isempty(w2)
   w2 = 1.1; 
end
if nargin < 4 || isempty(l)
   l = 4;
end
if nargin < 5 || isempty(m)
   m = 0.2; 
end
if nargin < 6 || isempty(depth)
   depth = 3;
end
% Checking if user-given inputs make sense. If not, set to default and
% throw a warning
if w1 >= w2 || w2 >= zmax || l >= 2*xmax || m <= 0 || l < 0 || w1 <= 0
   warning('Given inputs are incorrects. Set everything to default')
   w1 = 1;
   w2 = 1.1;
   l = 4;
   m = 0.2;
end

% Setting the right units of measure
w1 = w1*nm;
w2 = w2*nm;
l = l*nm;
depth = depth*eV;

% Preliminary steps
xstep = stepSize;
zstep = stepSize;

x_x = (-xmax:xstep:xmax)*nm;
z_z = (-zmax:zstep:zmax)*nm;

nx = length(x_x);
nz = length(z_z);

dx = min(diff(x_x));
dz = min(diff(z_z));

x = kron(x_x,z_z.^0);
z = kron(x_x.^0,z_z);

% potential U
U0 = depth;
U = zeros(size(x));
U(abs(z)<w1) = -U0;
U(abs(z)<w2 & abs(x)<l/2) = -U0;

% branch points
n = nx*nz;
Dxx_x = spdiags(ones(nx,1)*[1 -2 1]/dx^2,-1:1,nx,nx);
Dzz_z = spdiags(ones(nx,1)*[1 -2 1]/dz^2,-1:1,nz,nz);
H_L = -1/m*Dzz_z + diag(U(1:nz));
H_R = -1/m*Dzz_z + diag(U(end-nz+(1:nz)));

if norm(H_L - H_R) == 0
    % symmetric potential
    [V,D] = eig(H_L);
    [d,i] = sort(diag(D));
    V = V(:,i);
    p = zeros(nz,1); % 0 (left+right)
else
    % asymmetric potential
    [V_L,D_L] = eig(H_L);
    [V_R,D_R] = eig(H_R);
    V = [V_L,V_R];
    p = [-ones(nz,1);ones(nz,1)];
    [d,i] = sort([diag(D_L);diag(D_R)]);
    V = V(:,i);
    p = p(i); % -1 (left), 1 (right)
end

% matrix H
H = -1/m*(kron(Dxx_x,speye(nz))+kron(speye(nx),Dzz_z)) + spdiags(U',0,n,n);

% branch points
brpts = unique(d);

%% nonlinear matrices (SL = L_k, SU = U_k)
SL = cell(length(brpts),1);
if p(1) < 0
    SL{1} = [V(:,1);sparse(n-nz,1)];
elseif p(1) == 0
    SL{1} = [ [V(:,1);sparse(n-nz,1)], [sparse(n-nz,1);V(:,1)] ];
else
    SL{1} = [sparse(n-nz,1);V(:,1)];
end
c = 1;
for j = 2:length(d)
    if d(j-1) == d(j)
        if p(j) < 0
            SL{c} = [ SL{c}, [V(:,j);sparse(n-nz,1)] ];
        elseif p(j) == 0
            SL{c} = [ SL{c}, [V(:,j);sparse(n-nz,1)], [sparse(n-nz,1);V(:,j)] ];
        else
            SL{c} = [ SL{c}, [sparse(n-nz,1);V(:,j)] ];
        end
    else
        c = c + 1;
        if p(j) < 0
            SL{c} = [V(:,j);sparse(n-nz,1)];
        elseif p(j) == 0
            SL{c} = [ [V(:,j);sparse(n-nz,1)], [sparse(n-nz,1);V(:,j)] ];
        else
            SL{c} = [sparse(n-nz,1);V(:,j)];
        end
    end
end
SU = SL;
for j = 1:length(SL)
    SU{j} = SU{j}';
    SL{j} = -1/m/dx^2*SL{j};
end
% S_l = L_k*U_k
S = cell(length(brpts),1);
for j = 1:length(S)
    S{j} = SL{j}*SU{j};
end

%% nonlinear functions
f = cell(length(brpts),1);
for j = 1:interval-1
     f{j} = @(lambda) exp(1i.*(m.*(lambda-brpts(j))).^(0.5));
end
for j = interval:length(brpts)
     f{j} = @(lambda) exp(-(m.*(-lambda+brpts(j))).^(0.5));
end

% Setting coeffs
coeffs = cell(1, 2+nz);
coeffs{1} = H;
coeffs{2} = speye(n);
coeffs(3:end) = S(1:end);

% Setting xcoeffs
xcoeffs1 = SL';
xcoeffs2 = SU';
xcoeffs = {1 1 xcoeffs1{:}; H speye(n) xcoeffs2{:}};
% Setting fun and F
fun = @(lam) canyon_particle_fun(lam, f);
F = @(lam) canyon_particle_F(lam, xcoeffs, f);

end


function fun = canyon_particle_fun(lam, f)
lam = lam(:);
N = length(f);
fun = zeros(length(lam), N+2);
fun(:,1:2) = [ones(length(lam),1) -lam];
for k = 1:N
    fun(:,k+2) = -f{k}(lam);
end  
end

function F = canyon_particle_F(lam, xcoeffs, f)
N = length(f);
F = xcoeffs{2,1} - lam*xcoeffs{2,2};
for k = 1:N
    F = F  - f{k}(lam)*xcoeffs{1,k+2}*xcoeffs{2,k+2};
end

end

