% Seek [Ma0, R0, rho1, u1, p1] producing requested conditions at (R0, delta).
% Returns struct s with fields s.Ma0, s.R0, s.rho1, s.u1, and s.p1 producing a
% flow with observed values s.Ma_e, s.p_exi, and s.a2_e.
% Curried function results noz(Ly) and qoi(Ly) compute the flow on (R0, 0) to
% (R0, Ly) and quantities of interest at (R0, Ly), respectively.
%
% Parameter 'burn_in_iterations' (>0) fixes Ma0 for some number of iterations.
% Parameter 'opt' may supply additional nonlin_residmin options using optimset.
% However, options 'MaxIter' and 'fixed' will be ignored.
function [s, noz, qoi] = nozzle_baseflow(delta, gam0, Ma_e, p_exi, a2_e,    ...
                                         burn_in_iterations = 10,           ...
                                         opt = optimset('Algorithm',        ...
                                                        'lm_svd_feasible'))

  % Relative residuals of observations vs targets for [Ma0; R0; rho1; u1; p1]
  tgt = [Ma_e; p_exi; a2_e];
  f   = @(x) (obs_vector(delta, gam0, x(1), x(2), x(3), x(4), x(5)) - tgt)./tgt;

  % Establish bounds for [Ma0; R0; rho1; u1; p1], a guess, and constraint(s)
  opt = optimset(opt, 'lbound', [eps; eps; eps; -inf; eps],
                      'ubound', [inf; inf; inf;  inf; inf]);
  p = [Ma_e; 10*delta; 1; NaN; 1];  % Guess for Ma_e, R0, rho1, p1
  if Ma_e < 1;                      % Guess for u1 per Ma_e
    p(4) = max (-1/p(1), -realsqrt(2 / p(1)**2 / (gam0 - 1) + 1));
  else
    p(4) = mean([1/p(1); +realsqrt(2 / p(1)**2 / (gam0 - 1) + 1)]);
  end
  realizable = @(x) 2 / p(1)**2 / (gam0 - 1) + 1 - p(4)**2;
  opt = optimset(opt, 'inequc', {zeros(length(p)), ones(size(p)), realizable});

  % Solve the problem converting relative residual vector into absolute results
  % Fixes density and pressure and solves for other parameters in two phases
  pkg load odepkg optim;
  phase1 = optimset(opt, 'fixed', [1;0;1;0;1], 'MaxIter', burn_in_iterations);
  phase2 = optimset(opt, 'fixed', [0;0;1;0;1], 'MaxIter',               5000);
  try
    niter = 0; res = nan(3,1);
    [p, ans, cvg, outp] = nonlin_residmin(f, p, phase1); niter += outp.niter;
    [p, res, cvg, outp] = nonlin_residmin(f, p, phase2); niter += outp.niter;
  catch
    warning('nozzle_baseflow(%g, %g, %g, %g, %g) fails: %s', ...
            delta, gam0, Ma_e, p_exi, a2_e, lasterror.message);
    cvg = 0; p = nan(5,1);
  end
  res2 = norm(res.*optimget(phase2, 'weights', ones(size(tgt))), 2);
  res  = res.*tgt + tgt;

  % Package up problem specification, solution, results, and runtime behavior
  s = struct('delta', delta, 'gam0', gam0,
             'Ma0', p(1), 'R0', p(2), 'rho1', p(3), 'u1', p(4), 'p1', p(5),
             'Ma_e', res(1), 'p_exi', res(2), 'a2_e', res(3),
             'res2', res2, 'cvg', cvg, 'niter', niter);
  noz = @(Ly) nozzle    (s.Ma0, s.gam0, s.R0, realsqrt(s.R0**2+Ly**2), ...
                         s.u1, s.rho1, s.p1);
  qoi = @(Ly) nozzle_qoi(Ly, s.gam0, s.Ma0, s.R0, s.rho1, s.u1, s.p1);
end

% Repackage nozzle_qoi multiple return values into a vector of observations.
function f = obs_vector(delta, gam0, Ma0, R0, rho1, u1, p1)
  [Ma_e, p_exi, a2_e, a2_w] = nozzle_qoi(delta, gam0, Ma0, R0, rho1, u1, p1);
  f = [Ma_e; p_exi; a2_e];  % Ignore a2_w
end

%!test
%! opt=optimset('debug', 1);  % Case A
%! tic(), s=nozzle_baseflow(1, 1.4000, 0.4000, -0.0200, 1.2000, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));

%!test
%! opt=optimset('debug', 1);  % Case B
%! tic(), s=nozzle_baseflow(1, 1.4080, 0.9825, -0.0099, 4.2952, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));

%!test
%! opt=optimset('debug', 1);  % Case C
%! tic(), s=nozzle_baseflow(1, 1.4083, 1.1094, -0.0097, 4.3205, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));

%!test
%! opt=optimset('debug', 1);  % Case D
%! tic(), s=nozzle_baseflow(1, 1.4081, 1.0482, -0.0097, 4.3134, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));

%!test
%! opt=optimset('debug', 1);  % Case E
%! tic(), s=nozzle_baseflow(1, 1.4091, 0.4112, -0.0179, 4.1291, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));

%!test
%! opt=optimset('debug', 1);  % Case F
%! tic(), s=nozzle_baseflow(1, 1.4095, 0.1217, -0.0958, 3.9670, 10, opt), toc()
%! assert(s.res2 < sqrt(eps));
