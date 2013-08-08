% ODE integrator-ready implementation of velocity derivative u'(r, u)
up = @(r,u,Ma,g0) -(u./r) .* (2 + Ma**2*(g0-1) - Ma**2*(g0-1)*u.**2) ...
                          ./ (2 + Ma**2*(g0-1) - Ma**2*(g0+1)*u.**2);

% Computations for sound speed squared, density, and pressure given r, u
a2  = @(u,Ma,g0) ...
      1 + Ma**2 * (g0 - 1) * (1 - u.**2) / 2;
rho = @(r,u,a2,Ma,g0,rho1) ...
      rho1 * exp(-2*pi*Ma**2 * cumtrapz(r, r.*u.*up(r,u,Ma,g0) ./ a2));
p   = @(r,u,rho,Ma,g0,p1) ...
      p1 - 2*pi*Ma**2 * cumtrapz(r, r.*rho.*u.*up(r,u,Ma,g0));

% Establish tolerances, integration options, domain size, and parameters
pkg load odepkg;
tol = sqrt(eps);
vopt = odeset('RelTol',      tol, 'AbsTol',  tol, ...
              'InitialStep', tol, 'MaxStep', sqrt(sqrt(tol)));
R1 = 1; R2 = 2;
Ma = 1; g0 = 1.4;

% Solve initial value problem for supersonic nozzle
[sup_r, sup_u] = ode45(up, [R1 R2], 1/Ma + tol,    vopt, Ma, g0);

% Postprocess to obtain thermodynamic state and then plot
sup_a2  = a2 (sup_u,                 Ma, g0   );
sup_rho = rho(sup_r, sup_u, sup_a2,  Ma, g0, 1);
sup_p   = p  (sup_r, sup_u, sup_rho, Ma, g0, 1);
figure();
plot(sup_r, sup_u, 'r-', sub_r, sup_rho, 'b-', sup_r, sup_p, 'g-');
legend('velocity', 'density', 'pressure', 'location', 'northwest');
title('Supersonic nozzle: inflow -> outflow')
xlabel('radius');

% Solve initial value problem for subsonic nozzle
[sub_r, sub_u] = ode45(up, [R1 R2], -1/Ma/(1+tol), vopt, Ma, g0);

% Postprocess to obtain thermodynamic state and then plot REVERSING X
sub_a2  = a2 (sub_u,                 Ma, g0   );
sub_rho = rho(sub_r, sub_u, sub_a2,  Ma, g0, 1);
sub_p   = p  (sub_r, sub_u, sub_rho, Ma, g0, 1);
figure();
plot(sub_r, -1*sub_u, 'r-', sub_r, sub_rho, 'b-', sub_r, sub_p, 'g-');
legend('velocity', 'density', 'pressure', 'location', 'northeast');
title('Subsonic nozzle: inflow -> outflow')
xlabel('radius'); set(gca,'XDir','Reverse')
