% Solve nozzle IVP for [u; rho; p] on [R1, R2] given [u1; rho1; p1]
% via a "coupled" ODE-based approach.  Plot results when no values requested.
% Default value for p1 satisfies ideal gas equation of state given rho1, gam.
function [r u rho p a2 up rhop pp] = nozzle(Ma, gam, R1, R2, u1, rho1,
                                            p1 = rho1/gam
                                               *(1+(gam-1)/2*Ma**2*(1-u1**2)),
                                            rtol=sqrt(eps), atol=sqrt(eps))
  [Ma2 gam1] = deal(Ma**2, gam-1);
  assert(u1**2 < 2 / Ma2 / gam1 + 1,
         'Ma=%g, gam=%g, u1=%g imply a**2 <= 0', Ma, gam, u1);
  vopt  = odeset('RelTol',      rtol,       'AbsTol',  atol,
                 'InitialStep', sqrt(atol), 'MaxStep', sqrt(sqrt(atol)));
  [r x] = ode45(@nozzle_rhs, [R1 R2], [u1 rho1 p1], vopt, Ma2, gam1);
  [u  rho  p    ] = deal(x(:,1), x(:,2), x(:,3));
  [up rhop pp a2] = nozzle_details(r, u, rho, p, Ma2, gam1);
  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end
end

% ODEs [u; rho; p]' given r, x=[u; rho; p], Ma2=Ma**2, gamm1=gam-1
function f = nozzle_rhs(r, x, Ma2, gamm1)
  [up, rhop, pp] = nozzle_details(r, x(1), x(2), x(3), Ma2, gamm1);
  f              = [up; rhop; pp];
end

% Find pointwise solution details given state, Ma2=Ma**2, gamm1=gam-1
function [up, rhop, pp, a2] = nozzle_details(r, u, rho, p, Ma2, gamm1)
  u2   = u.**2;
  C    = (2./Ma2 + gamm1.*(1 - u2));
  up   = (u.*C) ./ (r.*(2*u2 - C));
  pp   = -Ma2.*rho.*u.*up;
  a2   = 1 + 0.5*Ma2.*gamm1.*(1 - u2);
  rhop = pp ./ a2;
end

%!test
%! % Does a solution satisfy steady governing equations in radial setting?
%! % A verification test, including derivatives, against governing equations.
%! % Ideal gas EOS will be approximately satisfied for "large enough" radii.
%! pkg load odepkg; Ma=1.5; gam=1.4; Rin=10; Rout=Rin+1/2; u1=-2/7; rho1=9/10;
%! [r u rho p a2 up rhop pp] = nozzle(Ma, gam, Rin, Rout, u1, rho1);
%! assert(zeros(size(r))', (u.*rho./r+rho.*up+u.*rhop)', 11*eps);  # Mass
%! assert(pp', (-Ma.**2.*rho.*u.*up)', 10*eps);                    # Momentum
%! assert(a2', (1 + Ma.**2.*(gam-1)./2.*(1-u.**2))', 10*eps);      # Energy
%! assert((rho.*a2)', (gam.*p)', 10*eps);                          # Ideal EOS

%!demo
%! % Solve subsonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma = 1; gam = 1.4; Rin = 1; Rout = 2;
%! nozzle(Ma, gam, Rout, Rin, -2/7, 1, 1);
%! title('Subsonic nozzle');
%! print('nozzle_subsonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();

%!demo
%! % Subsonic cases may (more robustly) have nearly sonic outflows prescribed
%! pkg load odepkg; Ma = 1; gam = 1.4; Rin = 1; Rout = 2;
%! nozzle(Ma, gam, Rin, Rout, -1/Ma+sqrt(eps), 1, 1);

%!demo
%! % Solve supersonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma = 1; gam = 1.4; Rin = 1; Rout = 2;
%! nozzle(Ma, gam, Rin, Rout, 1/Ma+sqrt(eps), 1, 1);
%! title('Supersonic nozzle');
%! print('nozzle_supersonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();
