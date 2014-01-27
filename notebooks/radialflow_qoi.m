% Solve radial flow for Ma_e and p_exi at Cartesian point (R1, delta).
function [Ma_e p_exi] = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1)
  try
    [r u rho p a2 up rhop pp] = radialflow(Ma, gam, R1,
                                           sqrt(R1.^2+delta.^2), u1, rho1, p1);
    Ma_e  = Ma * r(1) * abs(u(end)) ...
          / (r(end) * realsqrt(a2(end)));
    p_exi = sign(u(end)) * r(end) * delta * pp(end) ...
          / (Ma.^2 * r(1) * rho(end) * u(end).^2);
  catch
    warning('radialflow_qoi(%g, %g, %g, %g, %g, %g, %g) fails: %s',
            delta, gam, Ma, R1, u1, rho1, p1, lasterror.message);
    Ma_e = p_exi = NaN;
  end
end
