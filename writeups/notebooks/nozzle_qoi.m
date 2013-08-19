% Solve a minimal nozzle problem computing Ma_e, p_exi, T_e at (R0, delta).
function [Ma_e, p_exi, T_e] = nozzle_qoi(delta, gam0, Ma, R0, rho1, u1)

  R1   = max(eps, R0);           % Defensive as BFGS algorithm may try R0 == 0
  R2   = sqrt(R0**2 + delta**2); % Minimal domain with edge as final result
  Ma_e = p_exi = T_e = NaN;      % Produce NaNs if nozzle(...) ODE solver fails
  p1   = 0;                      % Absolute pressure is irrelevant to results
  try
    [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1);
    Ma_e  = Ma * r(end) * abs(u(end)) / (R2 * sqrt(a2(end)));
    p_exi = sign(u(end)) * R2 * pp(end) / (R0 * Ma * Ma);
    T_e   = a2(end);
  catch
    warning('nozzle_qoi(%g, %g, %g, %g, %g, %g) NaNs', ...
            delta, gam0, Ma, R0, rho1, u1);
  end

end
