function [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1)
% Solve nozzle initial value problem for u(r) on [R1, R2] given u1(R1)
% via an ODE integrator-ready implementation of velocity derivative u'(r, u).
% Afterward, compute sound speed squared, density, and pressure from r, u
% as well as derivatives of u and p.  Plots results when no values requested.

  Ma2   = Ma*Ma;
  up    = @(r,u) -(u./r) .* (2 + Ma2*(gam0 - 1) - Ma2*(gam0 - 1)*u.**2) ...
                         ./ (2 + Ma2*(gam0 - 1) - Ma2*(gam0 + 1)*u.**2);
  tol   = sqrt(eps);
  vopt  = odeset('RelTol',      tol, 'AbsTol',  tol, ...
                 'InitialStep', tol, 'MaxStep', sqrt(sqrt(tol)));
  [r,u] = ode45(up, [R1 R2], u1, vopt);                           % Full order
  up    = up(r, u);                                               % Shadow 'up'
  a2    = 1 + 0.5*Ma2*(gam0 - 1)*(1 - u.**2);                     % Full order
  rho   = rho1 * exp(-2*pi*Ma2 * cumtrapz(r, r.*u.*up ./ a2));    % Lower order
  pp    = - Ma2 * rho.*u.*up;                                     % Lower order
  p     = p1 + 2*pi * cumtrapz(r, r.*pp);                         % Lower order

  if 0 == nargout
    figure();
    plot(r, u, 'r-', r, rho, 'b-', r, p, 'g-');
    legend('Velocity', 'Density', 'Pressure', 'location', 'northwest');
    xlabel('Radius');
  end

end % function
