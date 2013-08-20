% Driver aspiring to target Ma_e, p_exi, and a2_e quantities at (R0, delta).
% Returns struct s with parameters s.Ma, s.R0, r.rho1, and s.u1 creating
% a flow producing observed values s.obs.Ma_e, s.obs.p_exi, and s.obs.a2_e.
% Parameter 'pin' may be used to provide an initial guess [Ma; R0; rho1; u1].
function s = nozzle_baseflow(delta, gam0, Ma_e, p_exi, a2_e,             ...
                             tol = sqrt(eps), maxiter = 100, debug = 0, ...
                             pin = [])
  if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg optim; end

  % Relative residuals of observations against targets for [Ma, R0, rho1, u1]
  target = [Ma_e; p_exi; a2_e];
  f = @(x) (baseflow_f(delta, gam0, x(1), x(2), x(3), x(4)) - target)./target;

  % Establish initial guess as well as lower and upper bounds.
  % Guess for u1 = x(end) should be consistent with sub- vs supersonic Ma_e
  if isempty(pin); pin =   [1;       1;       1;       sign(Ma_e-1)*Ma_e+eps]; end
  opt = optimset('lbound', [eps;     eps;     eps;     -realmax             ],
                 'ubound', [realmax; realmax; realmax; +realmax             ],
                 'TolFun', tol, 'MaxIter', maxiter, 'debug', debug);

  % Solve the problem and convert relative residuals into optimization results
  [p, res, cvg, outp] = nonlin_residmin(f, pin, opt);
  res2 = norm(res, 2); res .*= target; res += target;

  % Pack problem specification, solution, and runtime behavior into a struct
  s = struct('delta', delta, 'gam0', gam0,
             'Ma', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4), 'res2', res2,
             'req', struct('Ma_e', Ma_e,   'p_exi', p_exi,  'a2_e', a2_e   ),
             'obs', struct('Ma_e', res(1), 'p_exi', res(2), 'a2_e', res(3)),
             'cvg', cvg, 'tol', tol, 'maxiter', maxiter, 'niter', outp.niter);

  % Curry so s.nozzle(Ly) computes results on segment (R0,0) to (R0,Ly)
  % and so that s.qoi(Ly) computes quantities of interest on that segment
  s.nozzle = @(Ly) nozzle(s.Ma,s.gam0,s.R0,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,0);
  s.qoi    = @(Ly) nozzle_qoi(Ly, s.gam0, s.Ma, s.R0, s.rho1, s.u1);

end

% Repackage nozzle_qoi multiple return values into a column vector.
function f = baseflow_f(delta, gam0, Ma, R0, rho1, u1)
  [Ma_e, p_exi, a2_e] = nozzle_qoi(delta, gam0, Ma, R0, rho1, u1);
  f = [Ma_e; p_exi; a2_e];
end

%!demo
%! % Find [Ma, R0, rho1, u1] producing given supersonic conditions at (R0,1)
%! tic(), s = nozzle_baseflow(1, 1.4, 1.4, -0.02, 4.2), toc()
%! % Now s.nozzle(2) computes baseflow suitable for use from (R0,0) to (R0,2)

%!demo
%! % Find [Ma, R0, rho1, u1] producing given subsonic conditions at (R0,1)
%! tic(), s = nozzle_baseflow(1, 1.4, 0.4, -0.02, 4.2), toc()
%! % Now s.nozzle(2) computes baseflow suitable for use from (R0,0) to (R0,2)
