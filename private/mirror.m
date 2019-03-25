function [coeffs,fun,F] = mirror(state)
%MIRROR    Quartic PEP from calibration of cadioptric vision system.
%  [COEFFS,FUN,F] = nlevp('mirror',S) returns the coefficients of
%  a 9-by-9 quartic matrix polynomial
%  lambda^4*A4 + lambda^3*A3 + lambda^2*A2 + lambda*A1 + A0
%  arising in a homography-based method for calibrating a central
%  cadioptric vision system.
%  S (default: S = 0) is used in the function call RAND('twister',S) to
%  seed the random number generator (see HELP RAND); if S < 0 then
%  the newer syntax RNG(-S) is used.
%  Set S = 'nostate' to leave the random number generator unseeded.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2, A3, A4}.
%  F is the function handle lambda^4*A4 + lambda^3*A3 + lambda^2*A2 + lambda*A1 + A0
%  This problem has the properties pep, real, random.

%  Reference:
%  B. Zhang and Y. F. Li. A method for calibrating the central catadioptric
%  camera via homographic matrix.  In Proceedings of the 2008 IEEE
%  International Conference on Information and Automation, Zhangjiajie,
%  China, pages 972-977, 2008.

%  Based on function supplied by Beiwei Zhang.

if nargin < 1, state = []; end

if isempty(state)
  warning('NLEVP:random','Random numbers are now seeded by this function.')
  state = 0;
end

if ~strcmpi(state,'noseed')
   rng(state);
end

num = 5; % Number of points.
a = rand(3,1);
Rx = [1 0 0;0 cos(a(1)) sin(a(1));0 -sin(a(1)) cos(a(1))];
Ry = [cos(a(2)) 0 -sin(a(2));0 1 0;sin(a(2)) 0 cos(a(2))];
Rz = [cos(a(3)) sin(a(3)) 0;-sin(a(3)) cos(a(3)) 0;0 0 1];
R = Rx*Ry*Rz;
T = 15*rand(3,1);

a = 3; % Which means the panoramic mirror is Z = (X^2+Y^2-a^2)/(2a).

N = rand(3,1);
Xp = rand(1,num);
Yp = rand(1,num);
Zp = (ones(1,num)-N(1)*Xp-N(2)*Yp)/N(3);
P = [Xp;Yp;Zp];

pm  = fprojection(a,eye(3),zeros(3,1),P);
pmm = fprojection(a,R,T,P);

A0 = zeros(10,9); A1 = A0; A2 = A0; A3 = A0; A4 = A0; 
for i=1:num

    u = pm(1,i); v = pm(2,i); t = u*u+v*v;
    x = pmm(1,i); y = pmm(2,i); s = x*x+y*y;

    Res10 = [0,  0,  -s*t,  0,  0,  0,  0,  0,  0];
    Res20 = [ 0,  0,  0,  0,  0,  -s*t,  0,  0,  0];
    Res11 = [-2*s*u,  -2*s*v,  0,  0,  0,  0,  0,  0,  2*x*t];
    Res21 = [ 0,  0,  0,  -2*s*u,  -2*s*v,  0,  0,  0,  2*y*t];
    Res12 = [0,  0,  s+t,  0,  0,  0,  4*x*u,  4*x*v,  0];
    Res22 = [ 0,  0,  0,  0,  0,  s+t,  4*y*u,  4*y*v,  0];
    Res13 = [2*u,  2*v,  0,  0,  0,  0,  0,  0,  -2*x];
    Res23 = [ 0,  0,  0,  2*u,  2*v,  0,  0,  0,  -2*y];
    Res14 = [0,  0,  -1,  0,  0,  0,  0,  0,  0];
    Res24 = [ 0,  0,  0,  0,  0,  -1,  0,  0,  0];

    A0(2*i-1:2*i,:) = [Res10;Res20];
    A1(2*i-1:2*i,:) = [Res11;Res21];
    A2(2*i-1:2*i,:) = [Res12;Res22];
    A3(2*i-1:2*i,:) = [Res13;Res23];
    A4(2*i-1:2*i,:) = [Res14;Res24];

end

coeffs = {A0(1:9,:),A1(1:9,:),A2(1:9,:),A3(1:9,:),A4(1:9,:)};
fun = @(lam) nlevp_monomials(lam,4);
F = @(lam) coeffs{1} + lam.*(coeffs{2} + lam.*(coeffs{3}+lam*(coeffs{4} +lam*coeffs{5})));
end


% minima = min(e_mirror'-a);
% index = find(e_mirror'-a==minima); % Find the needed value
% mir_cal=e_mirror(index);
% HH = reshape(e_Homo(:, index), 3,3);
% e_Homo = HH';
% Esti_matrix=e_Homo/e_Homo(3,3) % This is the estimated matrix.
% true_matrix=R+T*N';
% true_matrix/true_matrix(3,3) % This is the true matrix (should be equal to estimated matrix)

% The panoramic mirror is Z=(X^2+Y^2-a^2)/(2a).
function pm  = fprojection(a,R,t,Ps)
num = size(Ps,2);
pm = zeros(3,num);
for i=1:num
    p = R*Ps(:,i)+t;
    x = p(1); y = p(2); z = p(3);
    langta = a*(z+sqrt(x^2+y^2+z^2))/(x^2+y^2);
    pm(:,i) = langta*p;
end
end
