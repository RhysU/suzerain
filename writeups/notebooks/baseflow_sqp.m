function s = baseflow_sqp(p_exi, dstar, gam0, Ma_e, T_e, maxiter, tol)
% Driver aspiring to target p_exi, Ma_e, and T_e quantities at (R0, dstar).
% A sample invocation in the spirit of 4m leeward of the stagnation point:
%   s = baseflow_sqp(-0.05, 1, 1.4, 1.1, 4.2)
%
% Establishes a sequential quadratic programming problem and solve with sqp:
%    minimize phi(x) subject to h(x) >= 0, lb <= x <= ub
% where x may be unpacked using
%    Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4);
% Solution and access to the sqp solver behavior are returned.
% See 'help sqp'.  If not supplied, maxiter and tol default to sqp defaults.
  if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end
  if nargin < 6; maxiter = 100;       end
  if nargin < 7; tol     = sqrt(eps); end

  % Establish initial guess, lower bounds, and upper bounds for each parameter
  x   =zeros(4, 1); lb   =zeros(size(x)); ub   =zeros(size(x));
  x(1)=1;           lb(1)= eps;           ub(1)=realmax;  % Ma   positive
  x(2)=1;           lb(2)= eps;           ub(2)=realmax;  % R0   positive
  x(3)=1;           lb(3)= eps;           ub(3)=realmax;  % rho1 positive
  x(4)=NaN;         lb(4)=-realmax;       ub(4)=realmax;  % u1   unleashd

  % Guess for u1 = x(end) must be consistent with sub- vs supersonic Ma_e
  x(end) = sign(Ma_e - 1) * Ma_e + eps;

  % Run sequential quadratic programming collecting results into a struct
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
end

function h = baseflow_h(dstar, gam0, Ma_e, p_exi, T_e, x)
% Specify h(x) >= 0 constraints using eps to accomplish h(x) > 0 where needed.
% For subsonic, -1/Ma < u1 < 0 achieves at most nearly sonic outflow.
% For supersonic, u1 > 1/Ma achieves at least nearly sonic inflow.
  Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4);
  h = merge(Ma_e < 1, [1/Ma-u1-eps; u1-eps], [u1-1/Ma-eps]);
end

function phi = baseflow_phi(dstar, gam0, Ma_e, p_exi, T_e, x)
% Compute phi(x) returning squared norm of mismatch in Ma_e, p_exi, T_e.
  Ma = x(1); R0 = x(2); rho1 = x(3); u1 = x(4); x
  R1 = R0; R2 = sqrt(R0**2 + dstar**2);
  [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, 1);
  phi = (Ma_e  - Ma*r(end)*abs(u(end)) / (R2*sqrt(a2(end))))**2 ...
      + (p_exi + R2*abs(pp(end)) / (R0*Ma*Ma)              )**2 ...
      + (T_e   - a2(end)                                   )**2
end
