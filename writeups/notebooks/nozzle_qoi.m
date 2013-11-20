% Solve minimal nozzle for Ma_e, p_exi, a2_e at (R0,delta) and a2_w at (R0,0).
function [Ma_e p_exi a2_e a2_w] = nozzle_qoi(delta, gam0, Ma0, R0, rho1, u1, p1)
  try
    [r, u, rho, p, a2, up, rhop, pp] ...
          = nozzle(Ma0, gam0, max(R0,eps), realsqrt(R0**2+delta**2), u1, rho1, p1);
    Ma_e  = Ma0*r(1)*abs(u(end)) / (r(end)*realsqrt(a2(end)));
    p_exi = sign(u(end))*r(end)*delta*pp(end) / (Ma0**2*r(1)*rho(end)*u(end)**2);
    a2_e  = a2(end);
    a2_w  = a2(1);
  catch
    warning('nozzle_qoi(%g, %g, %g, %g, %g, %g, %g) fails: %s', ...
            delta, gam0, Ma0, R0, rho1, u1, p1, lasterror.message);
    Ma_e = p_exi = a2_e = a2_w = NaN;
  end
end
