% Solve baseflow IVP for [u; rho; p] given [gam H Ma r1 r2 u1 rho1 p1]
% via a "coupled" ODE-based approach.  Plot results when no values requested.
function [r u rho p T up rhop pp] = baseflow(gam, H, Ma r1, r2, u1, rho1, p1,
                                             rtol=sqrt(eps), atol=sqrt(eps))
  vopt  = odeset('RelTol',      rtol,       'AbsTol',  atol,
                 'InitialStep', sqrt(atol), 'MaxStep', sqrt(sqrt(atol)));
  [gam1 twoH Ma2] = deal(gam-1, 2*H, Ma^2);
  [r x          ] = ode45(@rhs, [r1 r2], [u1 rho1 p1], vopt, gam1, twoH, Ma2);
  [u  rho  p    ] = deal(x(:,1), x(:,2), x(:,3));
  [up rhop pp T ] = details(r, u, rho, p, gam1, twoH, Ma2);

  if 0 == nargout
    figure();
    plot(r, u, '-', r, rho, '-', r, p, '-', r, T, '-');
    legend('Velocity', 'Density', 'Pressure', 'Temperature', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end
end

% ODEs [u; rho; p]' given r, x=[u; rho; p], gam1=gam-1, twoH=2*H, Ma2=Ma.^2
function f = rhs(r, x, gam1, twoH, Ma2)
  [up rhop pp] = details(r, x(1), x(2), x(3), gam1, twoH, Ma2);
  f            = [up; rhop; pp];
end

% Find pointwise solution details given state, gam1=gam-1, twoH=2*H, Ma2=Ma.^2
function [up, rhop, pp, T] = details(r, u, rho, p, gam1, twoH, Ma2)
  Ma2u2 = Ma2*u.^2;
  up    = (u ./ r) .* (gam1*twoH       - gam1*Ma2u2)
                   ./ ((gam1+2).*Ma2u2 - gam1*twoH );
  T     = gam1 * (twoH - Ma2u2) / 2;
  pp    = -Ma2 * rho .* u .* up;
  rhop  = pp ./ T;
end
