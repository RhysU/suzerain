% Solve a nozzle problem for Ma_e, p_exi, a2_e at (R0,delta) and a2_w at (R0,0).
function [Ma_e, p_exi, a2_e, a2_w] = nozzle_qoi(delta, gam0, Ma, R0, rho1, u1, p1)

  R1   = max(eps, R0);               % Defensive as BFGS may try R0 == 0
  R2   = sqrt(R0**2 + delta**2);     % Minimal domain with edge as final result
  Ma_e = p_exi = a2_e = a2_w = NaN;  % Produce NaNs if nozzle(...) fails
  try
    [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1);
    Ma_e  = Ma * R2 * abs(u(end)) / (R2 * sqrt(a2(end)));
    p_exi = sign(u(end)) * R2 * pp(end) / (R0 * Ma * Ma);
    a2_e  = a2(end);
    a2_w  = a2(1);
  catch
    warning('nozzle_qoi(%g, %g, %g, %g, %g, %g) NaNs', ...
            delta, gam0, Ma, R0, rho1, u1);
  end

end
