% ODE integrator-ready implementation of velocity derivative u'(r, u)
up = @(r,u,Ma,g0) -(u./r) .* (2 + Ma**2*(g0-1) - Ma**2*(g0-1)*u.**2) ...
                          ./ (2 + Ma**2*(g0-1) - Ma**2*(g0+1)*u.**2);

% Computations for sound speed squared, density, and pressure given r, u
a2    = @(u,Ma,g0) 1 + Ma**2 * (g0 - 1) * (1 - u.**2) / 2;
rho   = @(r,u,a2,Ma,g0,rho1) ...
        rho1 * exp(-2*pi*Ma**2 * cumtrapz(r, r.*u.*up(r,u,Ma,g0) ./ a2));
gradp = @(r,u,rho,Ma,g0) Ma**2 * rho.*u.*up(r,u,Ma,g0);

% Establish tolerances, integration options, domain size, and parameters
tol  = sqrt(eps);
pkg load odepkg;
vopt = odeset("RelTol", tol, "AbsTol", tol, "InitialStep", sqrt(tol));
R1 = 1; R2 = 4;
Ma = 1; g0 = 1.4;

% Subsonic converging nozzle with outflow R1 and inflow R2
[sub_r, sub_u] = ode45(up, [R1 R2], -1/Ma/(1+tol), vopt, Ma, g0);
sub_a2    = a2(sub_u, Ma, g0);
sub_rho   = rho(sub_r, sub_u, sub_a2, Ma, g0, 1);
sub_gradp = gradp(sub_r, sub_u, sub_rho, Ma, g0);

% Supersonic diverging nozzle with inflow R1 and outflow R2
[sup_r, sup_u] = ode45(up, [R1 R2], 1/Ma + tol,    vopt, Ma, g0);
sup_a2  = a2(sup_u, Ma, g0);
sup_rho = rho(sup_r, sup_u, sup_a2, Ma, g0, 1);
sup_gradp = gradp(sup_r, sup_u, sup_rho, Ma, g0);

% Plot both results
plot(sub_r, sub_u, '-', sup_r, sup_u, '-');
xlabel('r'); ylabel('u');
legend('Subsonic', 'Supersonic', 'location', 'southoutside');
