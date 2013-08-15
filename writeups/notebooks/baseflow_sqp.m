% Driver aspiring to target Ma_e, p_exi, and T_e quantities at (R0, dstar).
% The invocation
%   s = baseflow_sqp(1, 1.4, 1.1, -0.01, 4.2)
% produces parameters s.Ma, s.R0, s.rho1, and s.u1 with observed behaviors
% s.obs.Ma_e, s.obs.p_exi, and s.obs.T_e.  The SQP solver results appear in
% s.x, s.obj, s.info, s.iter, and s.nf.  See 'help sqp' for solver details.
% When not supplied, maxiter = 100 and tol = eps.
function s = baseflow_sqp(dstar, gam0, Ma_e, p_exi, T_e, maxiter, tol)
  if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end
  if nargin < 6; maxiter = 100; end
  if nargin < 7; tol     = eps; end

  % Establish initial guess, lower bounds, and upper bounds for parameters
  %    Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4);
  % Guess for u1 = x(end) must be consistent with sub- vs supersonic Ma_e
  x  = [1,       1,       1,       sign(Ma_e-1)*Ma_e+eps];
  lb = [eps,     eps,     eps,     -realmax             ];
  ub = [realmax, realmax, realmax, realmax              ];

  % Run sequential quadratic programming collecting results into a struct:
  % That is, minimize phi(x) subject to h(x) >= 0, lb <= x <= ub
  % Curries baseflow_h and baseflow_phi to expose them as h(x) and phi(x)
  s = struct('dstar', dstar, 'gam0',  gam0,
             'Ma_e',  Ma_e, 'p_exi', p_exi, 'T_e', T_e,
             'h',   @(x) baseflow_h  (dstar, gam0, Ma_e, p_exi, T_e, x),
             'phi', @(x) baseflow_phi(dstar, gam0, Ma_e, p_exi, T_e, x));
  [s.x, s.obj, s.info, s.iter, s.nf] ...
        = sqp(x, s.phi, [], s.h, lb, ub, maxiter, tol);
  s.Ma = s.x(1); s.R0 = s.x(2); s.rho1 = s.x(3); s.u1 = s.x(4);

  % Curry so that s.nozzle(Ly) provides data on (R0,0) to (R0,Ly)
  s.nozzle=@(Ly) nozzle(s.Ma,s.gam0,s.R0,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,1);

  % Finally, compute the observed edge quantities at dstar
  s.obs = struct(); [s.obs.phi,s.obs.Ma_e,s.obs.p_exi,s.obs.T_e] = s.phi(s.x);
end

% Specify h(x) >= 0 constraints using eps to accomplish h(x) > 0 where needed.
% For subsonic, -1/Ma < u1 < 0 achieves at most nearly sonic outflow.
% For supersonic, u1 > 1/Ma achieves at least nearly sonic inflow.
function h = baseflow_h(dstar, gam0, Ma_e, p_exi, T_e, x)
  Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4);
  h = merge(Ma_e < 1, [1/Ma-u1-eps; -u1-eps], [u1-1/Ma-eps]);
end

% Compute phi(x) returning squared norm of mismatch vs tMa_e, tp_exi, tT_e.
% Defends against R0 == eps as BFGS seemingly ignores lb <= x <= ub.
function [phi, Ma_e, p_exi, T_e] = baseflow_phi(dstar, gam0, tMa_e, ...
                                                tp_exi, tT_e, x)
  Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4);
  R1 = max(eps, R0); R2 = sqrt(R0**2 + dstar**2);
  try
    [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, 1);
    Ma_e  =   Ma*r(end)*abs(u(end)) / (R2*sqrt(a2(end)));
    p_exi = - R2*abs(pp(end)) / (R0*Ma*Ma);
    T_e   =   a2(end);
    phi = (tMa_e - Ma_e)**2 + (tp_exi - p_exi)**2  + (tT_e - T_e)**2
  catch
    phi = realmax; Ma_e = p_exi = T_e = NaN;
  end
end
