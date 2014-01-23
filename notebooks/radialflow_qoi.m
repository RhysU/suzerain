% Solve radial flow for Ma_e, p_exi and approximately ideal T_e at (R1, delta).
function [Ma_e p_exi T_e] = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1, theta=1)
  try
    [r u rho p a2 up rhop pp] = radialflow(Ma, gam, R1, sqrt(R1.^2+delta.^2), ...
                                           u1, rho1, p1, theta);
    Ma_e  = Ma*r(1)*abs(u(end)) / (r(end)*realsqrt(a2(end)));
    p_exi = sign(u(end))*r(end)*delta*pp(end) / (Ma.^2*r(1)*rho(end)*u(end).^2);
    T_e   = gam * Ma_e.^2 * p(end) / Ma.^2 / rho(end);
  catch
    warning('radialflow_qoi(%g, %g, %g, %g, %g, %g, %g, %g) fails: %s', ...
            delta, gam, Ma, R1, u1, rho1, p1, theta, lasterror.message);
    Ma_e = p_exi = T_e = NaN;
  end
end
