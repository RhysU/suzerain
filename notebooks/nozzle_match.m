% Produce a flow matching the given boundary layer edge conditions.
% If supplied, T_e specifies an edge temperature for the ideal gas EOS.
function [Ma R0 R uR rhoR pR] = nozzle_match(delta, gam, Ma_e, p_exi, T_e=0)
  R0   = max(roots([
          - abs(Ma_e.^2 - 1)*abs(p_exi),
            delta,
          - Ma_e.^2 * delta.^2 * abs(p_exi) * sign(Ma_e.^2 - 1)
         ]));
  Ma   = 1 / realsqrt(1/Ma_e.^2 + (gam - 1)*delta.^2/R0.^2/2);
  R    = realsqrt(R0.^2 + delta.^2);
  uR   = - R / R0 * sign(p_exi*(Ma_e.^2 - 1));
  rhoR = 1;
  pR   = merge(T_e == 0, 1, T_e*Ma.^2*rhoR/gam/Ma_e.^2);
  if (iscomplex(R0))
    warning('nozzle_match(%g, %g, %g, %g, %g) has complex root: %g', ...
            delta, gam, Ma_e, p_exi, T_e, R0);
  endif
  if T_e >= gam*Ma_e.^2*(1/2 + 1/Ma.^2/(gam - 1))
    warning('nozzle_match(%g, %g, %g, %g, %g) is nonrealizable', ...
            delta, gam, Ma_e, p_exi, T_e);
  endif
end

%!test % Round trip: Supersonic nozzle with hot edge
%! delta = 1; gam = 1.4087;
%! [Ma R0 R uR rhoR pR] = nozzle_match(delta, gam, 1.1906, -0.025439, 4.0040);
%! [r u rho p]          = nozzle(Ma, gam, R, R0, uR, rhoR, pR);         % R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));       % R0
%! [Ma_e p_exi T_e]     = nozzle_qoi(delta, gam, Ma, R1, u1, rho1, p1); % R0->R
%! assert([Ma_e p_exi T_e], [1.1906, -0.025439, 4.0040], -sqrt(eps));

%!test % Round trip: Subsonic nozzle with hot edge
%! delta = 1; gam = 1.4088;
%! [Ma R0 R uR rhoR pR] = nozzle_match(delta, gam, 0.54927, -0.014755, 4.1541);
%! [r u rho p]          = nozzle(Ma, gam, R, R0, uR, rhoR, pR);         % R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));       % R0
%! [Ma_e p_exi T_e]     = nozzle_qoi(delta, gam, Ma, R1, u1, rho1, p1); % R0->R
%! assert([Ma_e p_exi T_e], [0.54927, -0.014755, 4.1541], -sqrt(eps));

%!test % Round trip: Supersonic diffuser with cold edge with non-unit delta
%! delta = 0.5; gam = 1.4;
%! [Ma R0 R uR rhoR pR] = nozzle_match(delta, gam, 1.5, +0.02, 0.50);
%! [r u rho p]          = nozzle(Ma, gam, R, R0, uR, rhoR, pR);         % R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));       % R0
%! [Ma_e p_exi T_e]     = nozzle_qoi(delta, gam, Ma, R1, u1, rho1, p1); % R0->R
%! assert([Ma_e p_exi T_e], [1.5, +0.02, 0.50], -sqrt(eps));

%!test % Round trip: Subsonic diffuser without prescribed edge temperature
%! delta = 1; gam = 1.4;
%! [Ma R0 R uR rhoR pR] = nozzle_match(delta, gam, 0.5, +0.015);
%! [r u rho p]          = nozzle(Ma, gam, R, R0, uR, rhoR, pR);         % R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));       % R0
%! [Ma_e p_exi T_e]     = nozzle_qoi(delta, gam, Ma, R1, u1, rho1, p1); % R0->R
%! assert([Ma_e p_exi pR], [0.5, +0.015, 1.0], -sqrt(eps));
