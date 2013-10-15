% [eta, f, fp, fpp] = blasius_solver(etaf=8.8, eta0=0, f0, fp0, fpp0, tol)
% Solve Blasius boundary layer over eta [eta_0, eta_f] given initial conditions
%
% That is, advance the problem
%   f''' + f f'' / 2
% subject to initial conditions f(0) = f0, f'(0) = fp0, f''(0) = fpp0.
% The value tol is used for the relative and absolute tolerances.
%
% If not supplied, initial conditions are taken from results tabulated in
% "Highly Accurate Solutions of the Blasius and {Falkner-Skan} Boundary Layer
% Equations via Convergence Acceleration" by B.  D. Ganapol
% (http://arxiv.org/format/1006.3888v1).

function [eta, f, fp, fpp] = blasius_solver(etaf, eta0, f0, fp0, fpp0, tol)

  if nargin < 6; tol  = sqrt(sqrt(eps)); end
  if nargin < 5; fpp0 = 3.320573362E-01; end
  if nargin < 4; fp0  = 0.000000000E+00; end
  if nargin < 3; f0   = 0.000000000E+00; end
  if nargin < 2; eta0 = 0.0;             end
  if nargin < 1; etaf = 8.8;             end

  pkg load odepkg;
  vopt = odeset('RelTol', tol, 'AbsTol', tol,
                'InitialStep', tol, 'MaxStep', sqrt(tol),
                'OutputFcn', @odeplot);

  [eta, y] = ode45(@blasius_solver_f, [eta0 etaf], [f0;fp0;fpp0], vopt);
  f        = y(:,1);
  fp       = y(:,2);
  fpp      = y(:,3);

  if 0 == nargout;
    figure(); plot(eta, fp, '-'); xlabel('eta'); ylabel('fp');
  end

end

% Compute Blasius equation RHS [f'; f''; f'''] given [f; f'; f''].
function rhs = blasius_solver_f(eta, f)
  rhs = [ f(2); f(3); -f(1)*f(3)/2 ];
end
