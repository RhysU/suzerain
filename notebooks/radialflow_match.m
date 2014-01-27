% Produce a radial flow matching the given boundary layer edge conditions.
function [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, Ma_e, p_exi)
  % Track candidates from zero, one, or two positive real roots for R0
  R0   = roots([ - abs(Ma_e.^2 - 1)*abs(p_exi),
                   delta,
                 - Ma_e.^2 * delta.^2 * abs(p_exi) * sign(Ma_e.^2 - 1) ]);
  R0   = sort(R0(arrayfun(@isreal,R0) & real(R0) > 0), 'descend');
  Ma   = 1 ./ realsqrt(1/Ma_e.^2 + (gam - 1)*delta.^2./R0.^2/2);
  R    = realsqrt(R0.^2 + delta.^2);
  uR   = - R ./ R0 * sign(p_exi.*(Ma_e.^2 - 1));
  rhoR = ones(size(R0));
  pR   = rhoR .* Ma.^2 ./ gam ./ Ma_e.^2;

  % Thin candidates down to realizable solution with largest R0
  [ok iok] = max((pR ./ rhoR) < (1 ./ (gam - 1) + Ma.^2 / 2));
  if ok
    R0 = R0(iok); Ma   = Ma  (iok); R  = R (iok);
    uR = uR(iok); rhoR = rhoR(iok); pR = pR(iok);
  else
    warning('radialflow_match(%g, %g, %g, %g) has no realizable solution',
            delta, gam, Ma_e, p_exi);
    Ma = R0 = R = uR = rhoR = pR = NaN;
  end
end

%!test % Round trip: Supersonic nozzle
%! delta=1; gam=1.4087;
%! [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, 1.1906, -0.025439);
%! [r u rho p]          = radialflow(Ma, gam, R, R0, uR, rhoR, pR);         %R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));           %R0
%! [Ma_e p_exi]         = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1); %R0->R
%! assert([Ma_e p_exi], [1.1906, -0.025439], -sqrt(eps));

%!test % Round trip: Subsonic nozzle
%! delta=1; gam=1.4088;
%! [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, 0.54927, -0.014755);
%! [r u rho p]          = radialflow(Ma, gam, R, R0, uR, rhoR, pR);         %R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));           %R0
%! [Ma_e p_exi]         = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1); %R0->R
%! assert([Ma_e p_exi], [0.54927, -0.014755], -sqrt(eps));

%!test % Round trip: Supersonic diffuser with non-unit delta
%! delta=0.5; gam=1.4;
%! [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, 1.5, +0.02);
%! [r u rho p]          = radialflow(Ma, gam, R, R0, uR, rhoR, pR);         %R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));           %R0
%! [Ma_e p_exi]         = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1); %R0->R
%! assert([Ma_e p_exi], [1.5, +0.02], -sqrt(eps));

%!test % Round trip: Subsonic diffuser
%! delta=1; gam=1.4;
%! [Ma R0 R uR rhoR pR] = radialflow_match(delta, gam, 0.5, +0.015);
%! [r u rho p]          = radialflow(Ma, gam, R, R0, uR, rhoR, pR);         %R->R0
%! [R1 u1 rho1 p1]      = deal(r(end), u(end), rho(end), p(end));           %R0
%! [Ma_e p_exi]         = radialflow_qoi(delta, gam, Ma, R1, u1, rho1, p1); %R0->R
%! assert([Ma_e p_exi], [0.5, +0.015], -sqrt(eps));
