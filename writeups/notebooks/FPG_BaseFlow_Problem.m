% ODE integrator-ready implementation of u'(r, u)
up = @(r,u,Ma,gamma) -(u/r) * (2 + Ma**2*(gamma-1) - Ma**2*(gamma-1) * u*u) ...
                            / (2 + Ma**2*(gamma-1) - Ma**2*(gamma+1) * u*u);

% Establish tolerances, integration options, domain size, and parameters
tol  = sqrt(eps);
pkg load odepkg;
vopt = odeset("RelTol", tol, "AbsTol", tol, "InitialStep", tol);
R1 = 1; R2 = 4;
Ma = 1; gamma = 1.4;

% Subsonic converging nozzle with outflow R1 and inflow R2
[sub_r, sub_u] = ode45(up, [R1 R2], -1/Ma/(1+tol), vopt, Ma, gamma);

% Supersonic diverging nozzle with inflow R1 and outflow R2
[sup_r, sup_u] = ode45(up, [R1 R2], 1/Ma + tol,    vopt, Ma, gamma);

% Plot both results
plot(sub_r, sub_u, '-', sup_r, sup_u, '-');
xlabel('r'); ylabel('u');
legend('Subsonic', 'Supersonic', 'location', 'southoutside');
