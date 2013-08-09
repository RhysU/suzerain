% Script file loading any packages necessary when using Octave
1; if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end

% Solve nozzle initial value problem for u(r) on [R1, R2] given u1(R1)
% via an ODE integrator-ready implementation of velocity derivative u'(r, u).
% Afterward, compute sound speed squared, density, and pressure from r, u
function [r, u, a2, rho, p] = nozzle(Ma, gamma, R1, R2, u1, rho1, p1)

  Ma2    = Ma*Ma;
  up     = @(r,u) -(u./r) .* (2 + Ma2*(gamma - 1) - Ma2*(gamma - 1)*u.**2) ...
                          ./ (2 + Ma2*(gamma - 1) - Ma2*(gamma + 1)*u.**2);
  tol    = sqrt(eps);
  vopt   = odeset('RelTol',      tol, 'AbsTol',  tol, ...
                  'InitialStep', tol, 'MaxStep', sqrt(sqrt(tol)));
  [r, u] = ode45(up, [R1 R2], u1, vopt);
  up     = up(r, u); % Shadow
  a2     = 1 + 0.5*Ma2*(gamma - 1)*(1 - u.**2);
  rho    = rho1 * exp(-2*pi*Ma2 * cumtrapz(r, r.*u.*up ./ a2));
  p      = p1 - 2*pi*Ma2 * cumtrapz(r, r.*rho.*u.*up);

endfunction

% Define parameters for nozzle cases of interest
Ma = 1; gamma0 = 1.4; R1 = 1; R2 = 2;

% Solve supersonic nozzle and plot solution to file
[r, u, a2, rho, p] = nozzle(Ma, gamma0, R1, R2, 1/Ma + sqrt(eps), 1, 1);
fig = figure(); set(fig, "visible", "off");
plot(r, u, 'r-', r, rho, 'b-', r, p, 'g-');
legend('velocity', 'density', 'pressure', 'location', 'northwest');
title('Supersonic nozzle: inflow -> outflow');
xlabel('radius');
print('nozzle_supersonic.eps', '-depsc2');
close(fig);

% Solve subsonic nozzle and plot solution to file REVERSING X
[r, u, a2, rho, p] = nozzle(Ma, gamma0, R1, R2, -1/Ma/(1+sqrt(eps)), 1, 1);
fig = figure(); set(fig, "visible", "off");
plot(r, -1*u, 'r-', r, rho, 'b-', r, p, 'g-');
legend('-velocity', 'density', 'pressure', 'location', 'northeast');
title('Subsonic nozzle: inflow -> outflow');
xlabel('radius'); set(gca,'XDir','Reverse');
print('nozzle_subsonic.eps', '-depsc2');
close(fig);
