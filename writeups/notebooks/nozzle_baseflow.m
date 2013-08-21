% Find [Ma, R0, rho1, u1, p1] producing given conditions at (R0, delta).
% Returns struct s with s.Ma, s.R0, r.rho1, s.u1, and s.p1 producing a flow
% with observed values s.obs.Ma_e, s.obs.p_exi, and s.obs.a2_e.  Curried
% functions s.noz(Ly) and s.qoi(Ly) compute the flow on (R0, 0) to (R0, Ly)
% and quantities of interest at (R0, Ly), respectively.
%
% When not empty, 'pin' provides an initial guess [Ma; R0; rho1; u1; p1].
% Parameter 'opt' may supply options to nonlin_residmin via optimset.
% Useful options include TolFun, MaxIter, debug, and weights.
function s = nozzle_baseflow(delta, gam0, Ma_e, p_exi, a2_e, ...
                             pin = [], opt = optimset('MaxIter', 50))
  pkg load odepkg optim;

  % Residuals of observations against targets for [Ma, R0, rho1, u1, p1]
  % Relative (instead of absolute) residuals used otherwise Ma_e swamps p_exi!
  tgt = [Ma_e; p_exi; a2_e];
  f   = @(x) (baseflow_f(delta, gam0, x(1), x(2), x(3), x(4), x(5)) - tgt)./tgt;

  % Establish initial guess, feasible bounds, and fix an order one pressure
  % Guess for u1 = x(end) should be consistent with sub- vs supersonic Ma_e
  if isempty(pin); pin =        [  1;   1;   1;  sign(Ma_e-1)*Ma_e+eps;   1]; end
  opt = optimset(opt, 'lbound', [eps; eps; eps; -inf                  ; eps],
                      'ubound', [inf; inf; inf;  inf                  ; inf],
                      'fixed',  [  0;   0;   0;  0                    ;   1]);

  % Solve the problem and convert relative residuals into optimization results
  [p, res, cvg, outp] = nonlin_residmin(f, pin, opt);
  res2 = norm(res.*optimget(opt, 'weights', ones(size(tgt))), 2);
  res  = res.*tgt + tgt;

  % Pack problem specification, solution, and runtime behavior into a struct
  s = struct('delta', delta, 'gam0', gam0,
             'Ma', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4), 'p1', p(5),
             'req', struct('Ma_e', Ma_e,   'p_exi', p_exi,  'a2_e', a2_e   ),
             'obs', struct('Ma_e', res(1), 'p_exi', res(2), 'a2_e', res(3)),
             'res2', res2, 'cvg', cvg, 'niter', outp.niter);
  s.noz = @(Ly) nozzle(s.Ma,s.gam0,s.R0,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,s.p1);
  s.qoi = @(Ly) nozzle_qoi(Ly, s.gam0, s.Ma, s.R0, s.rho1, s.u1, s.p1);

end

% Repackage nozzle_qoi multiple return values into a column vector.
function f = baseflow_f(delta, gam0, Ma, R0, rho1, u1, p1)
  [Ma_e, p_exi, a2_e, a2_w] = nozzle_qoi(delta, gam0, Ma, R0, rho1, u1, p1);
  f = [Ma_e; p_exi; a2_e];  % Ignore a2_w
end

%!demo
%! tic(), s = nozzle_baseflow(1, 1.4, 1.4, -0.02, 1.2), toc()

%!demo
%! opt = optimset('TolFun', eps, 'MaxIter', 100, 'debug', 1);
%! tic(), s = nozzle_baseflow(1, 1.4, 0.4, -0.02, 1.2, [], opt), toc()
