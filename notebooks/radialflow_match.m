% Produce a flow matching the given boundary layer edge conditions.
% If supplied, T_e specifies an edge temperature for the ideal gas EOS.
function [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, Ma_e, p_exi,
                                                 T_e=0, theta=1)

  % Track candidates from zero, one, or two positive real roots for R0
  R0   = roots([ - abs(Ma_e.^2 - 1)*abs(p_exi),
                   delta,
                 - Ma_e.^2 * delta.^2 * abs(p_exi) * sign(Ma_e.^2 - 1) ]);
  R0   = sort(R0(arrayfun(@isreal,R0) & real(R0) > 0), 'descend');
  Ma   = 1 ./ realsqrt((1/Ma_e.^2 + (gam - 1)*delta.^2./R0.^2/2)./theta);
  R    = realsqrt(R0.^2 + delta.^2);
  uR   = - R ./ R0 * sign(p_exi.*(Ma_e.^2 - 1));
  rhoR = ones(size(R0));
  pR   = merge(T_e == 0, ones(size(R0)), T_e.*Ma.^2.*rhoR./gam./Ma_e.^2);

  % Thin candidates down to realizable solution with largest R0
  [ok iok] = max(T_e < gam.*Ma_e.^2.*(1/2 + theta./Ma.^2./(gam - 1)));
  if ok
    R0 = R0(iok); Ma   = Ma  (iok); R  = R (iok);
    uR = uR(iok); rhoR = rhoR(iok); pR = pR(iok);
  else
    warning('radialflow_match(%g, %g, %g, %g, %g) has no realizable solution',
            delta, gam, Ma_e, p_exi, T_e);
    Ma = R0 = R = uR = rhoR = pR = NaN;
  end

end

%!test % Round trip: Supersonic nozzle with hot edge
%! d=1; g=1.4087; t=1;
%! [Ma R0 R uR rhoR pR] = radialflow_match(d, g, 1.1906, -0.025439, 4.0040, t);
%! [r u rho p]          = radialflow(Ma, g, R, R0, uR, rhoR, pR, t);
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));
%! [Ma_e p_exi T_e]     = radialflow_qoi(d, g, Ma, R1, u1, rho1, p1, t);
%! assert([Ma_e p_exi T_e], [1.1906, -0.025439, 4.0040], -sqrt(eps));

%!test % Round trip: Subsonic nozzle with less hot edge
%! d=1; g=1.4088; t=1;
%! [Ma R0 R uR rhoR pR] = radialflow_match(d, g, 0.54927, -0.014755, 4.1541, t);
%! [r u rho p]          = radialflow(Ma, g, R, R0, uR, rhoR, pR, t);
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));
%! [Ma_e p_exi T_e]     = radialflow_qoi(d, g, Ma, R1, u1, rho1, p1, t);
%! assert([Ma_e p_exi T_e], [0.54927, -0.014755, 3.1541], -sqrt(eps));

%!test % Round trip: Supersonic diffuser with cold edge with non-unit delta
%! d=0.5; g=1.4; t=1;
%! [Ma R0 R uR rhoR pR] = radialflow_match(d, g, 1.5, +0.02, 0.50, t);
%! [r u rho p]          = radialflow(Ma, g, R, R0, uR, rhoR, pR, t);
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));
%! [Ma_e p_exi T_e]     = radialflow_qoi(d, g, Ma, R1, u1, rho1, p1, t);
%! assert([Ma_e p_exi T_e], [1.5, +0.02, 0.50], -sqrt(eps));

%!test % Round trip: Subsonic diffuser without prescribed edge temperature
%! d=1; g=1.4; t=1;
%! [Ma R0 R uR rhoR pR] = radialflow_match(d, g, 0.5, +0.015, 0, t);
%! [r u rho p]          = radialflow(Ma, g, R, R0, uR, rhoR, pR, t);
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));
%! [Ma_e p_exi T_e]     = radialflow_qoi(d, g, Ma, R1, u1, rho1, p1, t);
%! assert([Ma_e p_exi pR], [0.5, +0.015, 1.0], -sqrt(eps));
