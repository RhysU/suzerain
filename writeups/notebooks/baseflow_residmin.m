% Driver aspiring to target Ma_e, p_exi, and T_e quantities at (R0, dstar).
function s = baseflow_residmin(dstar, gam0, Ma_e, p_exi, T_e, ...
                               tol = sqrt(eps), maxiter = 100, debug = 0)
  if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg optim; end

  % Relative residuals of observations against targets for [Ma, R0, rho1, u1]
  target = [Ma_e; p_exi; T_e];
  f = @(x) (baseflow_f(dstar, gam0, x(1), x(2), x(3), x(4)) - target)./target;

  % Establish guess, lower/upper bounds, and nonlinear inequality constraints
  % Guess for u1 = x(end) must be consistent with sub- vs supersonic Ma_e
  pin =                    [1;       1;       1;       sign(Ma_e-1)*Ma_e+eps];
  opt = optimset('lbound', [eps;     eps;     eps;     -realmax             ],
                 'ubound', [realmax; realmax; realmax; +realmax             ],
                 'inequc', { zeros(4), ones(4,1), ... % NOP linear constraint
                             @(x) baseflow_h(Ma_e, x(1), x(2), x(3), x(4))  },
                 'TolFun', tol, 'MaxIter', maxiter, 'debug', debug);

  % Solve the problem and convert relative residuals into optimization results
  [p, res, cvg, outp] = nonlin_residmin(f, pin, opt);
  res .*= target; res += target;

  % Pack the problem specification, solution, and behavior into a struct
  s = struct('dstar', dstar, 'gam0', gam0,
             'Ma', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4),
             'req', struct('Ma_e', Ma_e,   'p_exi', p_exi,  'T_e', T_e   ),
             'obs', struct('Ma_e', res(1), 'p_exi', res(2), 'T_e', res(3)),
             'cvg', cvg, 'niter', outp.niter);

  % Curry so that s.nozzle(Ly) provides data on (R0,0) to (R0,Ly)
  s.nozzle=@(Ly) nozzle(s.Ma,s.gam0,s.R0,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,1);

end

% Specify h(x) >= 0 constraints using eps to accomplish h(x) > 0.
function h = baseflow_h(Ma_e, Ma, R0, rho1, u1)
  if Ma_e < 1
    h = [1/Ma-u1-eps; -u1-eps]; % -1/Ma < u1 < 0 achieves at most sonic outflow
  else
    h = [u1-1/Ma-eps];          % u1 > 1/Ma achieves at least sonic inflow
  end
end

% Compute vector of nozzle_edge results containing [tMa_e; tp_exi; tT_e].
function f = baseflow_f(dstar, gam0, Ma, R0, rho1, u1)
  [Ma_e, p_exi, T_e] = nozzle_edge(dstar, gam0, Ma, R0, rho1, u1);
  f = [Ma_e; p_exi; T_e];
end
