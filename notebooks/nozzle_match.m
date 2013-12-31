% Produce a flow matching the given boundary layer edge conditions.
% If supplied, T_e specifies an edge temperature for the ideal gas EOS.
function [R0 Ma0 u1 rho1 p1] = nozzle_match(delta, gam0, Ma_e, p_exi, T_e=0)
  R0   = max(roots([
          - abs(Ma_e.^2 - 1)*abs(p_exi),
            delta,
          - Ma_e.^2 * delta.^2 * abs(p_exi) * sign(Ma_e.^2 - 1)
         ]));
  Ma0  = 1 / realsqrt(1/Ma_e.^2 + (gam0 - 1)*delta.^2/R0.^2/2);
  R    = realsqrt(R0.^2 + delta.^2);
  u1   = - R / R0 * sign(p_exi*(Ma_e.^2 - 1));
  rho1 = 1;
  p1   = merge(T_e == 0, 1, T_e*Ma0.^2*rho1/gam0/Ma_e.^2);
  if (iscomplex(R0))
    warning('nozzle_match(%g, %g, %g, %g, %g) has complex root: %g', ...
            delta, gam0, Ma_e, p_exi, T_e, R0);
  endif
end

%!test % Supersonic nozzle with hot edge
%! [R0 Ma0 u1 rho1 p1] = nozzle_match(1, 1.4087, 1.1906, -0.025439, 4.0040);
%! [Ma_e p_exi T_e]    = nozzle_qoi  (1, 1.4087, Ma0, R0, u1, rho1, p1);
%! assert(Ma_e,   1.1906,   -0.002);
%! assert(p_exi, -0.025439, -0.002);
%! assert(T_e,    4.0040,   -0.002);

%!test % Subsonic nozzle with hot edge
%! [R0 Ma0 u1 rho1 p1] = nozzle_match(1, 1.4088, 0.54927, -0.014755, 4.1541);
%! [Ma_e p_exi T_e]    = nozzle_qoi  (1, 1.4088, Ma0, R0, u1, rho1, p1);
%! assert(Ma_e,   0.54927,  -0.001);
%! assert(p_exi, -0.014755, -0.001);
%! assert(T_e,    4.1541,   -0.001);

%!test % Supersonic diffuser with cold edge with non-unit delta
%! [R0 Ma0 u1 rho1 p1] = nozzle_match(0.5, 1.4, 1.5, +0.02, 0.50);
%! [Ma_e p_exi T_e]    = nozzle_qoi  (0.5, 1.4, Ma0, R0, u1, rho1, p1);
%! assert(Ma_e,   1.5,  -0.0015);
%! assert(p_exi, +0.02, -0.0015);
%! assert(T_e,    0.5,  -0.0015);

%!test % Subsonic diffuser without prescribed edge temperature
%! [R0 Ma0 u1 rho1 p1] = nozzle_match(1, 1.4, 0.5, +0.015);
%! [Ma_e p_exi T_e]    = nozzle_qoi  (1, 1.4, Ma0, R0, u1, rho1, p1);
%! assert(Ma_e,   0.5,   -0.0001);
%! assert(p_exi, +0.015, -0.0001);
