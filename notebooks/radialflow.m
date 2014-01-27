% Solve radialflow IVP for [u; rho; p] given [Ma gam R1 R2 u1 rho1 p1]
% via a "coupled" ODE-based approach.  Plot results when no values requested.
function [r u rho p a2 up rhop pp] = radialflow(Ma, gam, R1, R2, u1, rho1, p1,
                                                tol=sqrt(eps))

  [Ma2 gam1] = deal(Ma.^2, gam-1);
  assert(u1.^2 < 2 / Ma2 / gam1 + 1,
         'Ma=%g, gam=%g, u1=%g imply a.^2 <= 0', Ma, gam, u1);

  vopt  = odeset('RelTol',  tol, 'InitialStep', 0.01*abs(R1-R2),
                 'AbsTol',  eps, 'MaxStep',     0.10*abs(R1-R2));
  [r x] = ode45(@radialflow_rhs, [R1 R2], [u1 rho1 p1], vopt, Ma2, gam1);
  [u  rho  p    ] = deal(x(:,1), x(:,2), x(:,3));
  [up rhop pp a2] = radialflow_details(r, u, rho, p, Ma2, gam1);

  if 0 == nargout
    figure();
    plot(r, u, 'o-', r, rho, '+-', r, p, 'x-', r, Ma*abs(u)./sqrt(a2), '*-');
    legend('Velocity', 'Density', 'Pressure', 'Local Mach', ...
           'location', 'westoutside', 'orientation', 'vertical');
    xlabel('Radius');
    box('off');
  end
end

% ODEs [u; rho; p]' given r, x=[u; rho; p], Ma2=Ma.^2, gam1=gam-1
function f = radialflow_rhs(r, x, Ma2, gam1)
  [up, rhop, pp] = radialflow_details(r, x(1), x(2), x(3), Ma2, gam1);
  f              = [up; rhop; pp];
end

% Find pointwise derivatives and sound speed given state, Ma2=Ma.^2, gam1=gam-1
function [up, rhop, pp, a2] = radialflow_details(r, u, rho, p, Ma2, gam1)
  u2   = u.^2;
  C    = (2./Ma2 + gam1.*(1 - u2));
  up   = (u.*C) ./ (r.*(2*u2 - C));
  pp   = -Ma2.*rho.*u.*up;
  a2   = 1 + 0.5*Ma2.*gam1.*(1 - u2);
  rhop = pp ./ a2;
end

%!test
%! % Does a solution satisfy steady governing equations in radial setting?
%! % A verification test, including derivatives, against governing equations.
%! % Pressure p1 computed from ideal gas equation of state.
%! pkg load odepkg; Ma=1.5; gam=1.4; Rin=10; Rout=Rin+1/2; u1=-2/7; rho1=9/10;
%! p1 = rho1/gam *(1+(gam-1)/2*Ma.^2*(1-u1.^2));
%! [r u rho p a2 up rhop pp] = radialflow(Ma, gam, Rin, Rout, u1, rho1, p1);
%! assert(zeros(size(r))', (u.*rho./r+rho.*up+u.*rhop)', 10*eps);  # Mass
%! assert(pp', (-Ma.^2.*rho.*u.*up)', 10*eps);                     # Momentum
%! assert(a2', (1 + Ma.^2.*(gam-1)./2.*(1-u.^2))', 10*eps);        # Energy
%! assert((rho.*a2)', (gam.*p)', 10*eps);                          # Ideal EOS

%!demo % Solve subsonic nozzle (specifying outflow) and plot to file
%! pkg load odepkg; Ma=2; gam=1.4; Rin=1; Rout=Rin+1;
%! u_sonic = sqrt((2/Ma.^2 + gam - 1) / (gam + 1));
%! radialflow(Ma, gam, Rout, Rin, -u_sonic/5, 1, 1/2);
%! print('nozzle_subsonic.eps', '-depsc2', '-S640,480', '-F:9');
%! close();

%!demo % Solve supersonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma=1; gam=1.4; Rin=1; Rout=Rin+1;
%! u_sonic = sqrt((2/Ma.^2 + gam - 1) / (gam + 1));
%! radialflow(Ma, gam, Rin, Rout, 1.5*u_sonic, 1/2, 1);
%! print('nozzle_supersonic.eps', '-depsc2', '-S640,480', '-F:9');
%! close();
