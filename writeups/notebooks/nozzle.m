% Solve nozzle IVP for [u; rho; p] on [R1, R2] given [u1; rho1; p1]
% via a coupled ODE-based approach.  Plot results when no values requested.
function [r, u, rho, p, a2, up, pp] = nozzle(Ma0, gam0, R1, R2, u1, rho1, p1, ...
                                             reltol=sqrt(eps), abstol=sqrt(eps))

  assert(u1**2 < 2 / Ma0**2 / (gam0-1) + 1,
         'Ma0=%g, gam0=%g, u1=%g imply a**2 <= 0', Ma0, gam0, u1);
  vopt     = odeset('RelTol',      reltol,       'AbsTol',  abstol,
                    'InitialStep', sqrt(abstol), 'MaxStep', sqrt(sqrt(abstol)));
  [r, x]   = ode45(@nozzle_f, [R1 R2], [u1 log(rho1) p1], vopt, Ma0**2, gam0-1);
  u        =     x(:,1) ;
  rho      = exp(x(:,2));
  p        =     x(:,3) ;
  [up, a2] = nozzle_upa2(r, u, Ma0**2, gam0-1);
  pp       = - Ma0**2 * rho.*u.*up;

  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end

end

% Find [u; log rho; p]' given r, x=[u; log rho; p], Ma02=Ma0**2, gam0m1=gam0-1
function f = nozzle_f(r, x, Ma02, gam0m1)
  u = x(1); rho = exp(x(2)); p = x(3);        % Unpack
  [up, a2] = nozzle_upa2(r, u, Ma02, gam0m1); % Compute
  logrhop  = -Ma02*r*u*up / a2;
  pp       = -Ma02*r*rho*u*up;
  f = [up; logrhop; pp];                      % Pack
end

% Compute u' and a2 given r, u, Ma02=Ma0**2, gam0m1=gam0-1
function [up, a2] = nozzle_upa2(r, u, Ma02, gam0m1)
  up = -(u./r) .* (2 + Ma02*gam0m1 - Ma02*(gam0m1  )*u.**2) ...
               ./ (2 + Ma02*gam0m1 - Ma02*(gam0m1+2)*u.**2);
  a2 = 1 + 0.5*Ma02*gam0m1*(1 - u.**2);
end

%!demo
%! % Used identically to nozzle1 but produces more accurate solutions
%! pkg load odepkg; Ma0 = 1; gam0 = 1.4; Rinner = 1; Router = 2;
%! nozzle(Ma0, gam0, Rinner, Router, -1/Ma0+sqrt(eps), 1, 1);

%!test
%! % Subsonic code-to-code verification against simpler nozzle1 logic
%! pkg load odepkg; Ma0 = 2; gam0 = 1.4; Rinner = 1.5; Router = 2;
%! [r, u, rho, p, a2, up, pp] = nozzle1(Ma0, gam0, Router, Rinner, -2/7, 0.9, 1.1);
%! expected = [r(end), u(end), rho(end), p(end), a2(end), up(end), pp(end)];
%! [r, u, rho, p, a2, up, pp] = nozzle (Ma0, gam0, Router, Rinner, -2/7, 0.9, 1.1);
%! observed = [r(end), u(end), rho(end), p(end), a2(end), up(end), pp(end)];
%! assert(norm(expected - observed) / norm(expected) < sqrt(sqrt(eps)))
