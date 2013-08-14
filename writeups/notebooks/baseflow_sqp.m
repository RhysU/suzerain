function s = baseflow_sqp(p_exi, dstar, gam0, Ma_e, T_e, maxiter, tol)
% Driver aspiring to target p_exi, Ma_e, and T_e quantities at (R0, dstar).
% A sample invocation in the spirit of 4m leeward of the stagnation point:
%   s = baseflow_sqp(-0.05, 1, 1.4, 1.1, 4.2, 100, sqrt(eps))
%
% Establishes a sequential quadratic programming problem and solve with sqp:
%    minimize phi(x) subject to h(x) >= 0, lb <= x <= ub
% where x is the (alphabetically-ordered) vector which may be unpacked with
%    [p_exi, dstar, gam0, Ma, Ma_e, R0, rho1, T_e, u1] = num2cell(x){:};
% and lb, ub are used to hold some parameters fixed.  See 'help sqp'.
  if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end

  % Establish initial guess, lower bounds, and upper bounds for each parameter
  x   =zeros(11, 1); lb   =zeros(size(x)); ub   =zeros(size(x));
  x(1)=p_exi;        lb(1)=p_exi;          ub(1)=p_exi;    % p_exi  constant
  x(2)=dstar;        lb(2)=dstar;          ub(2)=dstar;    % dstar constant
  x(3)=gam0;         lb(3)=gam0;           ub(3)=gam0;     % gam0  constant
  x(4)=1;            lb(4)=eps;            ub(4)=realmax;  % Ma    positive
  x(5)=Ma_e;         lb(5)=Ma_e;           ub(5)=Ma_e;     % Ma_e  constant
  x(6)=1;            lb(6)=eps;            ub(6)=realmax;  % R0    positive
  x(7)=1;            lb(7)=eps;            ub(7)=realmax;  % rho1  positive
  x(8)=T_e;          lb(8)=T_e;            ub(8)=T_e;      % T_e   constant
  x(9)=NaN;          lb(9)=-realmax;       ub(9)=realmax;  % u1    unleashd

  % Guess for u1 = x(9) must be consistent with sub- vs supersonic Ma_e
  x(9) = sign(Ma_e - 1) * Ma_e + eps;

  % Run sequential quadratic programming collecting results into a struct
  s = struct('h', @(x) baseflow_h(x), 'phi', @(x) baseflow_phi(x));  % Expose
  [s.x, s.obj, s.info, s.iter, s.nf]                               ...
        = sqp(x, @baseflow_phi, [], @baseflow_h, lb, ub, maxiter, tol);
  [s.p_exi,s.dstar,s.gam0,s.Ma,s.Ma_e,s.R0,s.rho1,s.T_e,s.u1]       ...
        = num2cell(s.x){:};

  % Curry so that s.nozzle(Ly) provides data on (R0,0) to (R0,Ly)
  s.nozzle=@(Ly) nozzle(s.Ma,s.gam0,s.R0,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,1);
end

function h = baseflow_h(x)
% Specify h(x) >= 0 constraints using eps to accomplish h(x) > 0 where needed.
% For subsonic, -1/Ma < u1 < 0 achieves at most nearly sonic outflow.
% For supersonic, u1 > 1/Ma achieves at least nearly sonic inflow.
  [p_exi, dstar, gam0, Ma, Ma_e, R0, rho1, T_e, u1] = num2cell(x){:};
  h = merge(Ma_e < 1, [1/Ma-u1-eps; u1-eps], [u1-1/Ma-eps]);
end

function phi = baseflow_phi(x)
% Compute phi(x) on smallest domain returning mismatch in p_exi, Ma_e, T_e.
  [p_exi, dstar, gam0, Ma, Ma_e, R0, rho1, T_e, u1] = num2cell(x){:};
  R2 = sqrt(R0**2 + dstar**2);
  [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R0, R2, u1, rho1, 1);
  phi = (Ma_e  - Ma*r(end)*abs(u(end)) / (R2*sqrt(a2(end))))**2 ...
      + (p_exi + R2*abs(pp(end)) / (R0*Ma*Ma)              )**2 ...
      + (T_e   - a2(end)                                   )**2
end
