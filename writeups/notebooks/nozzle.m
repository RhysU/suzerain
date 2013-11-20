% Solve nozzle IVP for [u; rho; p] on [R1, R2] given [u1; rho1; p1]
% via a "coupled" ODE-based approach.  Plot results when no values requested.
% Default value for p1 satisfies ideal gas equation of state given rho1, gam.
function [r u rho p a2 up rhop pp] = nozzle(Ma, gam, R1, R2, u1, rho1,
                                            p1 = rho1/gam
                                               *(1+(gam-1)/2*Ma**2*(1-u1**2)),
                                            rtol=sqrt(eps), atol=sqrt(eps))
  assert(u1**2 < 2 / Ma**2 / (gam-1) + 1,
         'Ma=%g, gam=%g, u1=%g imply a**2 <= 0', Ma, gam, u1);

  vopt   = odeset('RelTol',      rtol,       'AbsTol',  atol,
                  'InitialStep', sqrt(atol), 'MaxStep', sqrt(sqrt(atol)));
  [r, x] = ode45(@nozzle_ode, [R1 R2], [u1 log(rho1) p1], vopt, Ma**2, gam-1);

  u                 =     x(:,1) ;
  rho               = exp(x(:,2));
  p                 =     x(:,3) ;
  [up, a2, logrhop] = nozzle_helper(r, u, Ma.**2, gam-1);
  rhop              = logrhop.*rho;
  pp                = -Ma.**2.*rho.*u.*up;

  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end
end

% ODEs [u; log rho; p]' given r, x=[u; log rho; p], Ma2=Ma**2, gamm1=gam-1
function f = nozzle_ode(r, x, Ma2, gamm1)
  [up, a2, logrhop] = nozzle_helper(r, x(1), Ma2, gamm1);
  pp                = -Ma2.*exp(x(2)).*x(1).*up;
  f                 = [up; logrhop; pp];
end

% Helper computing pointwise details given r, u, Ma2=Ma**2, gamm1=gam-1
function [up, a2, logrhop] = nozzle_helper(r, u, Ma2, gamm1)
  C       = (2./Ma2 + gamm1.*(1-u.**2));
  up      = (u .* C) ./ (r .* (2*u.**2 - C));
  a2      = 1 + 0.5*Ma2.*gamm1.*(1 - u.**2);
  logrhop = -Ma2.*u.*up ./ a2;
end

%!test
%! % FIXME Test with Ma != 1.0
%! % Does a solution satisfy steady governing equations in radial setting?
%! % A verification test, including derivatives, against governing equations.
%! % Ideal gas EOS will be approximately satisfied for "large enough" radii.
%! pkg load odepkg; Ma=1.0; gam=1.4; Rin=5; Rout=Rin+1; u1=-2/7; rho1=9/10;
%! [r u rho p a2 up rhop pp] = nozzle(Ma, gam, Rout, Rin, u1, rho1);
%! assert(zeros(size(r))', (u.*rho./r+rho.*up+u.*rhop)', 10*eps);  # Mass
%! assert(pp', (-rho.*u.*up)', 10*eps);                            # Momentum
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
