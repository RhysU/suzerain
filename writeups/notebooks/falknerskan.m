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
