% Solve radial nozzle for Ma_e, p_exi and approximately ideal T_e at (R0,delta).
function [Ma_e p_exi T_e] = nozzle_qoi(delta, gam0, Ma0, R0, u1, rho1=1, p1=1)
  try
    [r u rho p a2 up rhop pp] = nozzle(Ma0, gam0, max(R0,eps), ...
                                       realsqrt(R0.^2+delta.^2), u1, rho1, p1);
    Ma_e  = Ma0*r(1)*abs(u(end)) / (r(end)*realsqrt(a2(end)));
    p_exi = sign(u(end))*r(end)*delta*pp(end) / (Ma0.^2*r(1)*rho(end)*u(end).^2);
    T_e   = gam0 * Ma_e.^2 * p(end) / Ma0.^2 / rho(end);
  catch
    warning('nozzle_qoi(%g, %g, %g, %g, %g, %g, %g) fails: %s', ...
            delta, gam0, Ma0, R0, u1, rho1, p1, lasterror.message);
    Ma_e = p_exi = T_e = NaN;
  end
end
