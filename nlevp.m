function varargout = nlevp(name,varargin)
%NLEVP   Collection of nonlinear eigenvalue problems.
%  [COEFFS,FUN,OUT3,OUT4,...] = NLEVP(NAME,ARG1,ARG2,...)
%    generates the matrices and functions defining the problem specified by
%    NAME (a case insensitive string).
%    ARG1, ARG2,... are problem-specific input arguments.
%    All problems are of the form
%      T(lambda)*x = 0
%    where
%      T(lambda)= f0(lambda)*A0 + f1(lambda)*A1 + ... + fk(lambda)*Ak.
%    The matrices A0, A1, ..., Ak are returned in a cell array:
%    COEFFS = {A0,...,Ak}.
%    FUN is a function handle that can be used to evaluate the functions
%    f1(lambda),...,fk(lambda).  For a scalar lambda,
%    F = FUN(lambda) returns a row vector containing
%      F = [f1(lambda), f2(lambda), ..., fk(lambda)].
%    If lambda is a column vector, FUN(lambda) returns a row per element in
%    lambda.
%    [F,FP] = FUN(lambda) also returns the derivatives
%      FP = [f1'(lambda), f2'(lambda), ..., fk'(lambda)].
%    [F,FP,FPP,FPPP,...] = FUN(lambda) also returns higher derivatives.
%    OUT3, OUT4, ... are additional problem-specific output arguments.
%    See the list below for the available problems.
%
%  PROBLEMS = NLEVP('query','problems') or NLEVP QUERY PROBLEMS
%    returns a cell array containing the names of all problems
%    in the collection.
%  NLEVP('help','name') or NLEVP HELP NAME
%    gives additional information on problem NAME, including number and
%    meaning of input/output arguments.
%  NLEVP('query','name') or NLEVP QUERY NAME
%    returns a cell array containing the properties of the problem NAME.
%  PROPERTIES = NLEVP('query','properties') or NLEVP QUERY PROPERTIES
%    returns the properties used to classify problems in the collection.
%  NLEVP('query',property1,property2,...) or NLEVP QUERY PROPERTY1 ...
%    lists the names of all problems having all the specified properties.
%
%  [T,TP,TPP,...] = NLEVP('eval',NAME,LAMBDA,ARG1,ARG2,...)
%  evaluates the matrix function T and its derivatives TP, TPP,...
%  for problem NAME at the scalar LAMBDA.
%
%  NLEVP('version') or NLEVP VERSION
%    prints version, release date, and number of problems
%    of the installed NLEVP collection.
%  V = NLEVP('version')
%    returns a structure V containing version information.
%    V consists of the fields v.number, v.date, and v.problemcount.
%
%  Available problems:
%
%  acoustic_wave_1d   Acoustic wave problem in 1 dimension.
%  acoustic_wave_2d   Acoustic wave problem in 2 dimensions.
%  bcc_traffic        QEP from stability analysis of chain of non-identical cars.
%                     under bilateral cruise control.
%  bent_beam          6-by-6 NEP from a bent beam model.
%  bicycle            2-by-2 QEP from the Whipple bicycle model.
%  bilby              5-by-5 QEP from Bilby population model.
%  buckling_plate     3-by-3 NEP from a buckling plate model.
%  butterfly          Quartic matrix polynomial with T-even structure.
%  canyon_particle    NEP from the Schrödinger equation on a canyon-like shape.
%  cd_player          QEP from model of CD player.
%  circular_piston    Sparse QEP from model of circular piston.
%  clamped_beam_1d    NEP from 1D clamped beam model with delayed feedback control.
%  closed_loop        2-by-2 QEP associated with closed-loop control system.
%  concrete           Sparse QEP from model of a concrete structure.
%  damped_beam        QEP from simply supported beam damped in the  middle.
%  damped_gyro        QEP of a damped gyroscopic system.
%  deformed_consensus n-by-n QEP from multi-agent systems theory.
%  dirac              QEP from Dirac operator.
%  disk_brake100      100-by-100 QEP from a disk brake model.
%  disk_brake4669     4669-by-4669 QEP from a disk brake model.
%  distributed_delay1 3-by-3 NEP from distributed delay system.
%  elastic_deform     QEP from elastic deformation of anisotropic material.
%  fiber              NEP from fiber optic design.
%  foundation         Sparse QEP from model of machine foundations.
%  gen_hyper2         Hyperbolic QEP constructed from prescribed eigenpairs.
%  gen_tantipal2      T-anti-palindromic QEP with eigenvalues on the unit
%                     circle.
%  gen_tpal2          T-palindromic QEP with prescribed eigenvalues on the
%                     unit circle.
%  gun                NEP from model of a radio-frequency gun cavity.
%  hadeler            NEP due to Hadeler.
%  hospital           QEP from model of Los Angeles Hospital building.
%  intersection       10-by-10 QEP from intersection of three surfaces.
%  loaded_string      REP from finite element model of a loaded vibrating
%                     string.
%  metal_strip        QEP related to stability of electronic model of metal
%                     strip.
%  mirror             Quartic PEP from calibration of cadioptric vision system.
%  mobile_manipulator QEP from model of 2-dimensional 3-link mobile manipulator.
%  nep1               2-by-2 basic NEP example.
%  nep2               3-by-3 basic NEP example.
%  nep3               NEP with weighted norm coefficients.
%  neuron_dde         2-by-2 NEP from a neural-network DDE.
%  omnicam1           9-by-9 QEP from model of omnidirectional camera.
%  omnicam2           15-by-15 QEP from model of omnidirectional camera.
%  orr_sommerfeld     Quartic PEP arising from Orr-Sommerfeld equation.
%  pdde_stability     QEP from stability analysis of discretized PDDE.
%  pdde_symmetric     n-by-n NEP from a partial delay differential equation.
%  photonic_crystal   REP from dG-FEM of wave propagation in a periodic
%                     medium.
%  pillbox_cavity     170562-by-170562 NEP from a RF pillbox cavity.
%  pillbox_small      20-by-20 NEP from a RF pillbox cavity.
%  planar_waveguide   Quartic PEP from planar waveguide.
%  plasma_drift       Cubic PEP arising in Tokamak reactor design.
%  power_plant        8-by-8 QEP from simplified nuclear power plant problem.
%  qep1               3-by-3 QEP with known eigensystem.
%  qep2               3-by-3 QEP with known, nontrivial Jordan structure.
%  qep3               3-by-3 parametrized QEP with known eigensystem.
%  qep4               3-by-4 QEP with known, nontrivial Jordan structure.
%  qep5               3-by-3 nonregular QEP with known Smith form.
%  railtrack          QEP from study of vibration of rail tracks.
%  railtrack_rep      REP from study of vibration of rail tracks.
%  railtrack2         Palindromic QEP from model of rail tracks.
%  railtrack2_rep     REP from model of rail tracks.
%  relative_pose_5pt  Cubic PEP from relative pose problem in computer vision.
%  relative_pose_6pt  QEP from relative pose problem in computer vision.
%  sandwich_beam      NEP from model of a clamped sandwich beam.
%  schrodinger        QEP from Schrodinger operator.
%  schrodinger_abc    NEP from Schrodinger equation with absorbing boundary condition.
%  shaft              QEP from model of a shaft on bearing supports with a
%                     damper.
%  sign1              QEP from rank-1 perturbation of sign operator.
%  sign2              QEP from rank-1 perturbation of 2*sin(x) + sign(x)
%                     operator.
%  sleeper            QEP modelling a railtrack resting on sleepers.
%  speaker_box        QEP from finite element model of speaker box.
%  spring             QEP from finite element model of damped mass-spring
%                     system.
%  spring_dashpot     QEP from model of spring/dashpot configuration.
%  square_root        Square root of a skew-symmetric matrix.
%  surveillance       27-by-20 QEP from surveillance camera callibration.
%  time_delay         3-by-3 NEP from a time-delay system.
%  time_delay2        2-by-2 NEP from a time-delay system.
%  time_delay3        NEP with high-variance-norm coefficients.
%  utrecht1331        QEP  1331-by-1331 QEP with singular A1.
%  wing               3-by-3 QEP from analysis of oscillations of a wing in
%                     an airstream.
%  wiresaw1           Gyroscopic system from vibration analysis of wiresaw.
%  wiresaw2           QEP from vibration analysis of wiresaw with viscous
%                     damping.
%
%  Examples:
%  coeffs = nlevp('railtrack')
%    generates the matrices defining the railtrack problem.
%  nlevp('help','railtrack')
%    prints the help text of the railtrack problem.
%  nlevp('query','railtrack')
%    prints the properties of the railtrack problem.
%
%  For example code to solve all polynomial eigenvalue problems (PEPs)
%  in this collection of dimension < 500 with MATLAB's POLYEIG
%  see NLEVP_EXAMPLE.M.

%  Reference:
%  T. Betcke, N. J. Higham, V. Mehrmann, C. Schroeder, and F. Tisseur.
%  NLEVP: A Collection of Nonlinear Eigenvalue Problems,
%  ACM Transactions on Mathematical Software, 39(2), pp7:1-7:28, 2013.


% Check inputs
if nargin < 1, error('Not enough input arguments'); end
if ~ischar(name), error('NAME must be a string'); end

name = lower(name);

if strcmp(name,'query')
    if nargin == 1
       error('Not enough input arguments')
    end
    [varargout{1:nargout}] = nlevp_query(varargin{:});
    return
end

if strcmp('string',name)
    name = 'spring';
    warning('NLEVP:string_renamed','Problem string has been renamed spring.')
end

if strcmp('version',name)
    [varargout{1:nargout}] = nlevp_version(varargin{:});
    return
end

switch name
    case 'help'
        nlevp_home = which('nlevp');
        nlevp_home = strrep(nlevp_home, 'nlevp.m', '');
        if nargin < 2
           help nlevp
        else
           if ~nlevp_isoctave
               % some nlevp are shadowed by MATLAB functions. This fixes it
               if ispc % we need to see if we are on Windows or not
               help(sprintf('%sprivate\\%s', nlevp_home, varargin{1}))
               else
               help(sprintf('%sprivate/%s', nlevp_home, varargin{1}))
               end
           else
               % Uglier code necessary for Octave.
               eval(['help ', varargin{1}]);
           end
        end
    case 'eval'
        [varargout{1:max(nargout,1)}] = nlevp_eval(varargin{:});
    otherwise
        [varargout{1:nargout}] = feval(name,varargin{:});
end

end
