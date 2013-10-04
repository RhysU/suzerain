% [eta, f, fp, fpp] = falknerskan(beta0 = 1, beta = 0, alpha = Blasius, etaf=10)
% Solve Falkner-Skan IVP over eta=[0, eta_f] given beta0, beta, alpha.
%
% That is, advance the problem
%   f''' + beta0 f f'' + beta (1 - (f')^2) = 0
% subject to initial conditions f(0) = 0, f'(0) = 0, f''(0) = alpha
% where, when not specified, the parameters default to the Blasius
% solution described by Ganapol referencing coefficients by Boyd.
%
% Beyond http://en.wikipedia.org/wiki/Blasius_boundary_layer also see "Highly
% Accurate Solutions of the Blasius and {Falkner-Skan} Boundary Layer Equations
% via Convergence Acceleration" by B. D. Ganapol
% (http://arxiv.org/format/1006.3888v1) as well as "The Blasius Function in the
% Complex Plane" by John P. Boyd
% (http://dx.doi.org/10.1080/10586458.1999.10504626) for background.


function [eta, f, fp, fpp] = falknerskan(beta0, beta, alpha, etaf,
                                         reltol=sqrt(eps)/10,
                                         abstol=sqrt(eps)/10)

  if nargin < 4; etaf  = 10.0;                end
  if nargin < 3; alpha = 0.33205733621519630; end
  if nargin < 2; beta  = 0.0;                 end
  if nargin < 1; beta0 = 1.0;                 end

  % FIXME
  "WARNING: Function falknerskan is currently broken."

  pkg load odepkg;
  vopt = odeset('RelTol',      reltol,       'AbsTol',  abstol,
                'InitialStep', sqrt(abstol), 'MaxStep', sqrt(abstol));
  [eta, y] = ode45(@falknerskan_f, [0 etaf], [0; 0; alpha], vopt, beta0, beta);
  f        = y(:,1);
  fp       = y(:,2);
  fpp      = y(:,3);

  if 0 == nargout;
    figure(); plot(eta, fp, '-'); xlabel('eta'); ylabel('fp');
  end

end

% Compute Falkner-Skan RHS [f'; f''; f'''] given [f; f'; f''].
function rhs = falknerskan_f(eta, f, beta0, beta)
  rhs = [ f(2); f(3); -beta0*f(1)*f(3) - beta*(1-f(2).**2) ];
end

% Blasius profile values tabulated by Ganapol in Table 3b
%  0.0E+00  0.000000000E+00  0.000000000E+00  3.320573362E-01
%  2.0E-01  6.640999715E-03  6.640779210E-02  3.319838371E-01
%  4.0E-01  2.655988402E-02  1.327641608E-01  3.314698442E-01
%  6.0E-01  5.973463750E-02  1.989372524E-01  3.300791276E-01
%  8.0E-01  1.061082208E-01  2.647091387E-01  3.273892701E-01
%  1.0E+00  1.655717258E-01  3.297800312E-01  3.230071167E-01
%  1.2E+00  2.379487173E-01  3.937761044E-01  3.165891911E-01
%  1.4E+00  3.229815738E-01  4.562617647E-01  3.078653918E-01
%  1.6E+00  4.203207655E-01  5.167567844E-01  2.966634615E-01
%  1.8E+00  5.295180377E-01  5.747581439E-01  2.829310173E-01
%  2.0E+00  6.500243699E-01  6.297657365E-01  2.667515457E-01
%  2.2E+00  7.811933370E-01  6.813103772E-01  2.483509132E-01
%  2.4E+00  9.222901256E-01  7.289819351E-01  2.280917607E-01
%  2.6E+00  1.072505977E+00  7.724550211E-01  2.064546268E-01
%  2.8E+00  1.230977302E+00  8.115096232E-01  1.840065939E-01
%  3.0E+00  1.396808231E+00  8.460444437E-01  1.613603195E-01
%  3.2E+00  1.569094960E+00  8.760814552E-01  1.391280556E-01
%  3.4E+00  1.746950094E+00  9.017612214E-01  1.178762461E-01
%  3.6E+00  1.929525170E+00  9.233296659E-01  9.808627878E-02
%  3.8E+00  2.116029817E+00  9.411179967E-01  8.012591814E-02
%  4.0E+00  2.305746418E+00  9.555182298E-01  6.423412109E-02
%  4.2E+00  2.498039663E+00  9.669570738E-01  5.051974749E-02
%  4.4E+00  2.692360938E+00  9.758708321E-01  3.897261085E-02
%  4.6E+00  2.888247990E+00  9.826835008E-01  2.948377201E-02
%  4.8E+00  3.085320655E+00  9.877895262E-01  2.187118635E-02
%  5.0E+00  3.283273665E+00  9.915419002E-01  1.590679869E-02
%  5.2E+00  3.481867612E+00  9.942455354E-01  1.134178897E-02
%  5.4E+00  3.680919063E+00  9.961553040E-01  7.927659815E-03
%  5.6E+00  3.880290678E+00  9.974777682E-01  5.431957680E-03
%  5.8E+00  4.079881939E+00  9.983754937E-01  3.648413667E-03
%  6.0E+00  4.279620923E+00  9.989728724E-01  2.402039844E-03
%  6.2E+00  4.479457297E+00  9.993625417E-01  1.550170691E-03
%  6.4E+00  4.679356615E+00  9.996117017E-01  9.806151170E-04
%  6.6E+00  4.879295811E+00  9.997678702E-01  6.080442648E-04
%  6.8E+00  5.079259772E+00  9.998638190E-01  3.695625701E-04
%  7.0E+00  5.279238811E+00  9.999216041E-01  2.201689553E-04
%  7.2E+00  5.479226847E+00  9.999557173E-01  1.285698072E-04
%  7.4E+00  5.679220147E+00  9.999754577E-01  7.359298339E-05
%  7.6E+00  5.879216466E+00  9.999866551E-01  4.129031111E-05
%  7.8E+00  6.079214481E+00  9.999928812E-01  2.270775140E-05
%  8.0E+00  6.279213431E+00  9.999962745E-01  1.224092624E-05
%  8.2E+00  6.479212887E+00  9.999980875E-01  6.467978611E-06
%  8.4E+00  6.679212609E+00  9.999990369E-01  3.349939753E-06
%  8.6E+00  6.879212471E+00  9.999995242E-01  1.700667989E-06
%  8.8E+00  7.079212403E+00  9.999997695E-01  8.462841214E-07
