% Solve nozzle IVP for [u; rho; p] on [R1, R2] given [u1; rho1; p1]
% via a "coupled" ODE-based approach.  Plot results when no values requested.
function [r, u, rho, p, a2, up, rhop, pp] ...
    = nozzle(Ma0, gam0, R1, R2, u1, rho1, p1, rtol=sqrt(eps), atol=sqrt(eps))

  assert(u1**2 < 2 / Ma0**2 / (gam0-1) + 1,
         'Ma0=%g, gam0=%g, u1=%g imply a**2 <= 0', Ma0, gam0, u1);

  vopt   = odeset('RelTol',      rtol,       'AbsTol',  atol,
                  'InitialStep', sqrt(atol), 'MaxStep', sqrt(sqrt(atol)));
  [r, x] = ode45(@nozzle_ode, [R1 R2], [u1 log(rho1) p1], vopt, Ma0**2, gam0-1);

  u                 =     x(:,1) ;
  rho               = exp(x(:,2));
  p                 =     x(:,3) ;
  [up, a2, logrhop] = nozzle_helper(r, u, Ma0.**2, gam0-1);
  rhop              = logrhop.*rho;
  pp                = -Ma0.**2.*rho.*u.*up;

  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end

end

% ODEs [u; log rho; p]' given r, x=[u; log rho; p], Ma02=Ma0**2, gam0m1=gam0-1
function f = nozzle_ode(r, x, Ma02, gam0m1)
  [up, a2, logrhop] = nozzle_helper(r, x(1), Ma02, gam0m1);
  pp                = -Ma02.*exp(x(2)).*x(1).*up;
  f                 = [up; logrhop; pp];
end

% Helper computing pointwise details given r, u, Ma02=Ma0**2, gam0m1=gam0-1
function [up, a2, logrhop] = nozzle_helper(r, u, Ma02, gam0m1)
  C       = (2./Ma02 + gam0m1.*(1-u.**2));
  up      = (u .* C) ./ (r .* (2*u.**2 - C));
  a2      = 1 + 0.5*Ma02.*gam0m1.*(1 - u.**2);
  logrhop = -Ma02.*u.*up ./ a2;
end

%!test
%! % Does a solution satisfy steady governing equations in radial setting?
%! % A verification test, including derivatives, against governing equations.
%! pkg load odepkg; Ma0 = 1; gam0 = 1.4; Rin = 1; Rout = 2;
%! [r u rho p a2 up rhop pp] = nozzle(Ma0, gam0, Rout, Rin, -2/7, 0.9, 1.1);
%! assert(zeros(size(r)), u.*rho./r + rho.*up + u.*rhop, 10*eps);  # Mass
%! assert(pp, -rho.*u.*up, 10*eps);                                # Momentum
%! assert(a2, 1 + Ma0.**2.*(gam0-1)./2.*(1-u.**2), 10*eps);        # Energy

%!demo
%! % Solve subsonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma0 = 1; gam0 = 1.4; Rin = 1; Rout = 2;
%! nozzle(Ma0, gam0, Rout, Rin, -2/7, 1, 1);
%! title('Subsonic nozzle');
%! print('nozzle_subsonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();

%!demo
%! % Subsonic cases may (more robustly) have nearly sonic outflows prescribed
%! pkg load odepkg; Ma0 = 1; gam0 = 1.4; Rin = 1; Rout = 2;
%! nozzle1(Ma0, gam0, Rin, Rout, -1/Ma0+sqrt(eps), 1, 1);

%!demo
%! % Solve supersonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma0 = 1; gam0 = 1.4; Rin = 1; Rout = 2;
%! nozzle(Ma0, gam0, Rin, Rout, 1/Ma0+sqrt(eps), 1, 1);
%! title('Supersonic nozzle');
%! print('nozzle_supersonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();
