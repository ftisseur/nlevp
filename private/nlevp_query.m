function names = nlevp_query(varargin)
%QUERY         Query helper function for NLEVP.
%  See HELP NLEVP for syntax.

% All properties recognized by the collection.
properties={...
    'rep'
    'pep'
    'qep'
    'nep'
    'nonregular'
    'nonsquare'
    'real'
    'symmetric'
    'hermitian'
    'hyperbolic'
    'elliptic'
    'overdamped'
    'proportionally-damped'
    't-palindromic'
    '*-palindromic'
    't-anti-palindromic'
    '*-anti-palindromic'
    't-even'
    '*-even'
    't-odd'
    '*-odd'
    'gyroscopic'
    'random'
    'sparse'
    'parameter-dependent'
    'scalable'
    'solution'
    'tridiagonal'
    'banded'
    'low-rank'
    };

if nargin == 1 && strcmpi(varargin{1},'properties')
   names = properties; return
end

% COINCIDENCE is an array with (i,j) element 1 if problem i has property j
% and 0 otherwise.

coincidence = zeros(1,length(properties));  % Will add rows below.
problems = [];

vararg1 = varargin{1};
quick_return = false;

% Set up problem properties.
% Put the problems in alphabetical order by name.

% NOTE: Structure X is used only for Octave compatibility,
% since nested functions don't inherit variables in Octave.
% For MATLAB can delete all instances of X.
names = {}; X.QR = quick_return; X.N = names; X.Pb = problems;
X.Pp = properties; X.C = coincidence; X.V = vararg1;

X = add_problem(X,'acoustic_wave_1d','pep','qep','symmetric','*-even',...
                  'parameter-dependent','scalable','sparse',...
                  'tridiagonal', 'banded','low-rank');
X = add_problem(X,'acoustic_wave_2d','pep','qep','symmetric','*-even',...
                  'parameter-dependent','scalable','sparse', 'banded',...
                  'low-rank');
X = add_problem(X,'bcc_traffic','pep','qep','real','parameter_dependent',...
                  'scalable','sparse','tridiagonal','banded');
X = add_problem(X,'bent_beam','nep','real');
X = add_problem(X,'bicycle','pep','qep','real','parameter-dependent');
X = add_problem(X,'bilby','pep','qep','real','parameter-dependent');
X = add_problem(X,'buckling_plate','nep','real','symmetric');
X = add_problem(X,'butterfly','pep','real','parameter-dependent','T-even',...
                  'scalable','sparse', 'banded'); %(8,8) diagonals
X = add_problem(X,'canyon_particle','nep','parameter-dependent',...
                  'scalable','sparse', 'banded', 'low-rank'); %(80,80) diags)
X = add_problem(X,'cd_player','pep','qep','real');
X = add_problem(X, 'clamped_beam_1d', 'nep','real', 'parameter-dependent', 'scalable', ...
                 'sparse', 'tridiagonal', 'banded', 'low-rank');
X = add_problem(X,'circular_piston','pep','qep','real','sparse');
X = add_problem(X,'closed_loop','pep','qep','real','parameter-dependent');
X = add_problem(X,'concrete','pep','qep','symmetric','parameter-dependent',...
                  'sparse', 'low-rank');
X = add_problem(X,'damped_beam','pep','qep','real','symmetric','scalable',...
                  'sparse', 'banded', 'low-rank'); %(3,3) diags PROBLEM DESCRIPTION
X = add_problem(X,'damped_gyro','pep', 'qep', 'real','parameter-dependent',...
                  'scalable','sparse', 'banded'); %(7,7)
X = add_problem(X,'deformed_consensus','pep','qep','real','symmetric',...
                  'parameter-dependent','scalable');
X = add_problem(X,'dirac','pep','qep','real','symmetric',...
                  'parameter-dependent','scalable');
X = add_problem(X,'disk_brake100','pep','qep','real','parameter-dependent');
X = add_problem(X,'disk_brake4669','pep','qep','real','sparse',...
                  'parameter-dependent');
X = add_problem(X,'distributed_delay1', 'nep', 'real', 'solution');
X = add_problem(X,'elastic_deform','pep','qep','real','T-even',...
                  'parameter-dependent','sparse','scalable','banded');
X = add_problem(X,'fiber','nep','sparse','solution', 'tridiagonal',...
                  'banded', 'low-rank');
X = add_problem(X,'foundation','pep','qep','symmetric','sparse');
X = add_problem(X,'gen_hyper2','pep','qep','real','symmetric','hyperbolic',...
                  'parameter-dependent','scalable','solution','random');
X = add_problem(X,'gen_tantipal2','pep','qep','real','T-anti-palindromic',...
                  'parameter-dependent', 'scalable', 'random');
X = add_problem(X,'gen_tpal2','pep','qep','real','T-palindromic',...
                  'parameter-dependent', 'scalable', 'random');
X = add_problem(X,'gun','nep','sparse', 'banded', 'low-rank'); %(843,843)
X = add_problem(X,'hadeler','nep','real','symmetric','scalable');
X = add_problem(X,'hospital','pep','qep','real');
X = add_problem(X,'intersection','pep','qep','real');
X = add_problem(X,'loaded_string','rep','real','symmetric','scalable',...
                  'parameter-dependent','sparse','tridiagonal','banded',...
                  'low-rank');            
X = add_problem(X,'metal_strip','pep','qep','real');
X = add_problem(X,'mirror','pep','real','random');
X = add_problem(X,'mobile_manipulator','pep','qep','real');
X = add_problem(X,'nep1','nep');
X = add_problem(X,'nep2','nep');
X = add_problem(X,'nep3','nep','scalable','parameter-dependent','random');
X = add_problem(X,'neuron_dde','nep','real','parameter-dependent');
X = add_problem(X,'omnicam1','pep','qep','real');
X = add_problem(X,'omnicam2','pep','qep','real');
X = add_problem(X,'orr_sommerfeld','pep','parameter-dependent', ...
                  'scalable');
X = add_problem(X,'pdde_stability','pep','qep','scalable','parameter-dependent',...
                  'sparse','symmetric', 'banded'); %(15,15)
X = add_problem(X,'pdde_symmetric','nep','real','symmetric','scalable',...
                  'sparse', 'banded');
X = add_problem(X,'photonic_crystal','rep','symmetric','sparse',...
                  'parameter-dependent','scalable', 'low-rank');
X = add_problem(X,'pillbox_cavity','nep','real','symmetric','sparse',...
                  'banded','low-rank'); %(31603)
X = add_problem(X,'pillbox_small','nep','symmetric');
X = add_problem(X,'planar_waveguide','pep','real','symmetric',...
                  'tridiagonal', 'banded', 'low-rank');
X = add_problem(X,'plasma_drift','pep');
X = add_problem(X,'power_plant','pep','qep','symmetric','parameter-dependent');
X = add_problem(X,'qep1','pep','qep','real','solution');
X = add_problem(X,'qep2','pep','qep','real','solution');
X = add_problem(X,'qep3','pep','qep','real','parameter-dependent','solution');
X = add_problem(X,'qep4','pep','qep','nonregular','nonsquare','real',...
                  'solution');
X = add_problem(X,'qep5','pep','qep','nonregular','real');
X = add_problem(X,'railtrack','pep','qep','t-palindromic','sparse',...
                  'low-rank');
X = add_problem(X,'railtrack_rep','rep','sparse', 'low-rank');
X = add_problem(X,'railtrack2','pep','qep','t-palindromic','sparse',...
                  'scalable','parameter-dependent', 'low-rank');
X = add_problem(X,'railtrack2_rep','rep','sparse',...
                  'scalable','parameter-dependent', 'low-rank');
X = add_problem(X,'relative_pose_5pt','pep','real');
X = add_problem(X,'relative_pose_6pt','pep','qep','real');
X = add_problem(X,'sandwich_beam','nep','sparse', 'banded'); %(6,6)
X = add_problem(X,'schrodinger','pep','qep','real','symmetric','sparse');
X = add_problem(X,'schrodinger_abc','nep','real','sparse', 'banded', ...
                  'low-rank','scalable'); 
X = add_problem(X,'shaft','pep','qep','real','symmetric','sparse',...
                  'banded', 'low-rank'); % (3,3)
X = add_problem(X,'sign1','pep','qep','hermitian','parameter-dependent',...
                  'scalable');
X = add_problem(X,'sign2','pep','qep','hermitian','parameter-dependent',...
                  'scalable');
X = add_problem(X,'sleeper','pep','qep','real','scalable','sparse',...
                  'symmetric','proportionally-damped','solution');
X = add_problem(X,'speaker_box','pep','qep','real','symmetric');
X = add_problem(X,'spring','pep','qep','real','symmetric',...
                  'proportionally-damped','parameter-dependent',...
                  'scalable','sparse', 'tridiagonal', 'banded');
X = add_problem(X,'spring_dashpot','pep','qep','real','parameter-dependent',...
                  'scalable','sparse','random', 'low-rank');
X = add_problem(X,'square_root','nep','real','sparse');
X = add_problem(X,'surveillance','pep','qep','real','nonsquare','nonregular');
X = add_problem(X,'time_delay','nep', 'real');
X = add_problem(X,'time_delay2','nep', 'real', 'parameter-dependent');
X = add_problem(X,'time_delay3','nep','real','scalable','parameter-dependent','random');
X = add_problem(X,'utrecht1331','pep','qep', 'sparse', 'banded'); %(132,132)
X = add_problem(X,'wing','pep','qep','real');
X = add_problem(X,'wiresaw1','pep','qep','real','t-even','gyroscopic',...
                  'parameter-dependent','scalable');
X = add_problem(X,'wiresaw2','pep','qep','real','parameter-dependent',...
                  'scalable');

quick_return = X.QR; names = X.N; problems = X.Pb;
properties = X.Pp; coincidence = X.C;

if quick_return == true, return, end

if strcmpi(vararg1,'problems')
   names = problems'; return
end

columns = zeros(1,0);
for i=1:nargin
    if ~ischar(varargin{i}), error('Properties must be strings'); end
    foundproperty = false;
    for j=1:length(properties)
        if strcmpi(varargin{i},properties{j})
            foundproperty = true;
            columns = [columns,j];
            break;
        end
    end
    if ~foundproperty
        warning('NLEVP:unknown_property','Unknown property');
    end
end
rows = all(coincidence(:,columns),2);

if nargout == 0 && sum(rows) == 0
    disp('There is no problem in the collection matching your properties.');
end

[names{1:sum(rows)}] = problems{rows};

       % Nested function
       function X = add_problem(X, name,varargin)
       %ADD_PROBLEM   Adds new problem to data structures.

       quick_return = X.QR; names = X.N; problems = X.Pb;
       properties = X.Pp; coincidence = X.C; vararg1 = X.V;

       % NJH: Quick return business is slightly inelegant, but from
       % within ADD_PROBLEM can only return to QUERY, not NLEVP.  Also
       % didn't want to add test after each call to ADD_PROBLEM in QUERY.

       if quick_return == true, return, end
       if strcmpi(name,vararg1)
          % Quick return if properties for a matrix requested.
          names = varargin';
          quick_return = true;
       end

       n = length(problems) + 1;
       problems{n} = name;
       match = 0;
       for i = 1:length(varargin)
          prop = varargin{i};
          for j = 1:length(properties)
              if strcmpi(prop, properties{j})
                 match = j;
                 break
              end
          end
          if match == 0
             error(['Internal error: Property ' prop ' is illegal'])
          end
          coincidence(n,match) = 1;
       end

       X.QR = quick_return; X.N = names; X.Pb = problems;
       X.Pp = properties; X.C = coincidence;

       end % of function

names = names(:);
end
