% Seek [Ma, R0, rho1, u1, p1] producing requested conditions at (R0, delta).
% Returns struct s with fields s.Ma, s.R0, s.rho1, s.u1, and s.p1 producing a
% flow with observed values s.obs.Ma_e, s.obs.p_exi, and s.obs.a2_e.  Curried
% functions s.noz(Ly) and s.qoi(Ly) compute the flow on (R0, 0) to (R0, Ly) and
% quantities of interest at (R0, Ly), respectively.
%
% Parameter 'opt' may supply additional nonlin_residmin options using optimset.
function s = nozzle_baseflow(delta, gam0, Ma_e, p_exi, a2_e, ...
                             opt = optimset('TolFun', eps, 'MaxIter', 100))

  % Relative residuals of observations vs targets for [Ma; R0; rho1; u1; p1]
  tgt = [Ma_e; p_exi; a2_e];
  f   = @(x) (obs_vector(delta, gam0, x(1), x(2), x(3), x(4), x(5)) - tgt)./tgt;

  % Establish bounds for [Ma; R0; rho1; u1; p1] and a reasonable initial guess
  opt = optimset(opt, 'lbound', [eps; eps; eps; -inf; eps],
                      'ubound', [inf; inf; inf;  inf; inf]);
  p = ones(5,1);                             % Guess 1 for R0, rho1, p1
  p(1) = realsqrt(p(2)**2+delta**2) / p(2);  % Feasible guess for Ma per a2_e
  if Ma_e < 1;                               % Feasible guess for u1 per Ma_e
    p(4) = max (-1/p(1), -realsqrt(2 / p(1)**2 / (gam0 - 1) + 1));  % FIXME
  else
    p(4) = mean([1/p(1); -realsqrt(2 / p(1)**2 / (gam0 - 1) + 1)]); % FIXME
  end

  % Constrain x(1) = Ma and x(2) = R0 per to permit desired a2_e behavior
  if a2_e <= 1
    inequc = @(x)  [Ma_e * realsqrt(x(2)**2+delta**2) / x(2) - x(1)];
  else
    inequc = @(x) -[Ma_e * realsqrt(x(2)**2+delta**2) / x(2) - x(1)];
  end
  opt = optimset(opt, 'inequc', { zeros(length(p)), ones(size(p)), inequc });

  % Solve the problem converting relative residual vector into absolute results
  % Fixes density and pressure and solves for the remainder in phases
  pkg load odepkg optim;
  phase(1) = optimset(opt, 'fixed',  [ 1; 0; 1; 0; 1]);
  phase(2) = optimset(opt, 'fixed',  [ 0; 0; 1; 0; 1]);
  [p, res, cvg, outp] = nonlin_residmin(f, p, phase(1)); niter  = outp.niter;
  [p, res, cvg, outp] = nonlin_residmin(f, p, phase(2)); niter += outp.niter;
  res2 = norm(res.*optimget(phase(end), 'weights', ones(size(tgt))), 2);
  res  = res.*tgt + tgt;

  % Package up problem specification, solution, results, and runtime behavior
  s = struct('delta', delta, 'gam0', gam0,
             'Ma', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4), 'p1', p(5),
             'tgt', struct('Ma_e', Ma_e,   'p_exi', p_exi,  'a2_e', a2_e  ),
             'obs', struct('Ma_e', res(1), 'p_exi', res(2), 'a2_e', res(3)),
             'res2', res2, 'cvg', cvg, 'niter', niter);
  s.noz = @(Ly) nozzle    (s.Ma, s.gam0, s.R0, realsqrt(s.R0**2+Ly**2), ...
                           s.u1, s.rho1, s.p1);
  s.qoi = @(Ly) nozzle_qoi(Ly, s.gam0, s.Ma, s.R0, s.rho1, s.u1, s.p1);
end

% Repackage nozzle_qoi multiple return values into a vector of observations.
function f = obs_vector(delta, gam0, Ma, R0, rho1, u1, p1)
  [Ma_e, p_exi, a2_e, a2_w] = nozzle_qoi(delta, gam0, Ma, R0, rho1, u1, p1);
  f = [Ma_e; p_exi; a2_e];  % Ignore a2_w
end

%!demo
%! opt = optimset('TolFun', eps, 'MaxIter', 100, 'debug', 1);
%! tic(), s = nozzle_baseflow(1, 1.4, 0.4, -0.02, 1.2, opt), toc()

%!test
%! s=nozzle_baseflow(1, 1.4080, 0.9825, -0.0099, 4.2952), assert(s.res2 < 1e-2);

%!test
%! s=nozzle_baseflow(1, 1.4083, 1.1094, -0.0097, 4.3205), assert(s.res2 < 1e-2);

%!test
%! s=nozzle_baseflow(1, 1.4081, 1.0482, -0.0097, 4.3134), assert(s.res2 < 1e-2);

%!test
%! s=nozzle_baseflow(1, 1.4091, 0.4112, -0.0179, 4.1291), assert(s.res2 < 1e-2);

%!test
%! s=nozzle_baseflow(1, 1.4095, 0.1217, -0.0958, 3.9670), assert(s.res2 < 1e-2);
