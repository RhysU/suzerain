% Solve nozzle initial value problem for u(r) on [R1, R2] given u1(R1)
% via a scalar ODE-based approach.  Plot results when no values requested.
function [r, u, rho, p, a2, up, pp] = nozzle1(Ma, gam0, R1, R2, u1, rho1, p1, ...
                                              reltol=sqrt(eps), abstol=sqrt(eps))

  Ma2    = Ma*Ma;
  gam0m1 = gam0 - 1;
  up     = @(r,u) -(u./r) .* (2 + Ma2*gam0m1 - Ma2*(gam0m1  )*u.**2) ...
                          ./ (2 + Ma2*gam0m1 - Ma2*(gam0m1+2)*u.**2);
  vopt   = odeset('RelTol',      reltol,       'AbsTol',  abstol,
                  'InitialStep', sqrt(abstol), 'MaxStep', sqrt(sqrt(abstol)));
  [r,u]  = ode45(up, [R1 R2], u1, vopt);                    % Full accuracy
  a2     = 1 + 0.5*Ma2*gam0m1*(1 - u.**2);                  % Full accuracy
  up     = up(r, u);                                        % Shadow prior 'up'
  rho    = rho1 * exp(-Ma2 * cumtrapz(r, r.*u.*up ./ a2));  % Lower accuracy
  pp     = - Ma2 * rho.*u.*up;                              % Lower accuracy
  p      = p1 + cumtrapz(r, r.*pp);                         % Lower accuracy

  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end

end

%!demo
%! % Solve subsonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma = 1; gam0 = 1.4; Rinner = 1; Router = 2;
%! nozzle1(Ma, gam0, Router, Rinner, -2/7, 1, 1);  % Breaks for u1 = -1/3
%! title('Subsonic nozzle');
%! print('nozzle_subsonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();

%!demo
%! % Subsonic cases may (more robustly) have nearly sonic outflows prescribed
%! pkg load odepkg; Ma = 1; gam0 = 1.4; Rinner = 1; Router = 2;
%! nozzle1(Ma, gam0, Rinner, Router, -1/Ma+sqrt(eps), 1, 1);

%!demo
%! % Solve supersonic nozzle (specifying inflow) and plot to file
%! pkg load odepkg; Ma = 1; gam0 = 1.4; Rinner = 1; Router = 2;
%! nozzle1(Ma, gam0, Rinner, Router,  1/Ma+sqrt(eps), 1, 1);
%! title('Supersonic nozzle');
%! print('nozzle_supersonic.eps', '-depsc2', '-S512,384', '-F:8');
%! close();
