% Driver aspiring to target Ma_e, p_exi, and a2_e quantities at (R0, delta).
% Returns struct s with parameters s.Ma, s.R0, r.rho1, and s.u1 creating
% a flow producing observed values s.obs.Ma_e, s.obs.p_exi, and s.obs.a2_e.
%
% Parameter 'opt' supplies options to nonlin_residmin per optimset.
% Useful options include TolFun, MaxIter, debug, and weights.
% When not empty, 'pin' provides an initial guess [Ma; R0; rho1; u1].
function s = nozzle_baseflow(delta, gam0, Ma_e, p_exi, a2_e, ...
                             pin = [], opt = optimset('MaxIter', 50))
  pkg load odepkg optim;

  % Residuals of observations against targets for [Ma, R0, rho1, u1]
  % Relative (instead of absolute) residuals used otherwise Ma_e swamps p_exi!
  target = [Ma_e; p_exi; a2_e];
  f = @(x) (baseflow_f(delta, gam0, x(1), x(2), x(3), x(4)) - target)./target;

  % Establish initial guess as well as lower and upper bounds.
  % Guess for u1 = x(end) should be consistent with sub- vs supersonic Ma_e
  if isempty(pin); pin =        [1;   1;   1;    sign(Ma_e-1)*Ma_e+eps]; end
  opt = optimset(opt, 'lbound', [eps; eps; eps; -inf                  ],
                      'ubound', [inf; inf; inf;  inf                  ]);

  % Solve the problem and convert relative residuals into optimization results
  [p, res, cvg, outp] = nonlin_residmin(f, pin, opt);
  res2 = norm(res.*optimget(opt, 'weights', ones(size(target))), 2);
  res  = res.*target + target;

  % Pack problem specification, solution, and runtime behavior into a struct
  s = struct('delta', delta, 'gam0', gam0,
             'Ma', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4),
             'req', struct('Ma_e', Ma_e,   'p_exi', p_exi,  'a2_e', a2_e   ),
             'obs', struct('Ma_e', res(1), 'p_exi', res(2), 'a2_e', res(3)),
             'res2', res2, 'cvg', cvg, 'niter', outp.niter);

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
%! % Afterwards, s.nozzle(2) computes baseflow for use from (R0,0) to (R0,2)
%! tic(), s = nozzle_baseflow(1, 1.4, 1.4, -0.02, 4.2), toc()

%!demo
%! % Find [Ma, R0, rho1, u1] producing given subsonic conditions at (R0,1)
%! opt = optimset('TolFun', eps, 'MaxIter', 100, 'debug', 1);
%! tic(), s = nozzle_baseflow(1, 1.4, 0.4, -0.02, 4.2, [], opt), toc()
