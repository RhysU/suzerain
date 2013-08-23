% Solve nozzle IVP for [u; rho; p] on [R1, R2] given [u1; rho1; p1]
% via a coupled ODE-based approach.  Plot results when no values requested.
function [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1, ...
                                             reltol=sqrt(eps), abstol=sqrt(eps))

  assert(u1**2 < 2 / Ma**2 / (gam0-1) + 1,
         'Ma=%g, gam0=%g, u1=%g imply a**2 <= 0', Ma, gam0, u1);
  vopt     = odeset('RelTol',      reltol,       'AbsTol',  abstol,
                    'InitialStep', sqrt(abstol), 'MaxStep', sqrt(sqrt(abstol)));
  [r, x]   = ode45(@nozzle_f, [R1 R2], [u1 log(rho1) p1], vopt, Ma**2, gam0-1);
  u        =     x(:,1) ;
  rho      = exp(x(:,2));
  p        =     x(:,3) ;
  [up, a2] = nozzle_upa2(r, u, Ma**2, gam0-1);
  pp       = - Ma**2 * rho.*u.*up;

  if 0 == nargout
    figure();
    plot(r, u, '-', r, sqrt(a2), '-', r, rho, '-', r, p, 'g-');
    legend('Velocity', 'Sound speed', 'Density', 'Pressure', ...
           'location', 'north', 'orientation', 'horizontal');
    xlabel('Radius');
  end

end

% Find [u; log rho; p]' given r, x=[u; log rho; p], Ma2=Ma**2, gam0m1=gam0-1
function f = nozzle_f(r, x, Ma2, gam0m1)
  u = x(1); rho = exp(x(2)); p = x(3);       % Unpack
  [up, a2] = nozzle_upa2(r, u, Ma2, gam0m1); % Compute
  logrhop  = -Ma2*r*u*up / a2;
  pp       = -Ma2*r*rho*u*up;
  f = [up; logrhop; pp];                     % Pack
end

% Compute u' and a2 given r, u, Ma2=Ma**2, gam0m1=gam0-1
function [up, a2] = nozzle_upa2(r, u, Ma2, gam0m1)
  up = -(u./r) .* (2 + Ma2*gam0m1 - Ma2*(gam0m1  )*u.**2) ...
               ./ (2 + Ma2*gam0m1 - Ma2*(gam0m1+2)*u.**2);
  a2 = 1 + 0.5*Ma2*gam0m1*(1 - u.**2);
end

%!demo
%! % Used identically to nozzle1 but produces more accurate solutions
%! pkg load odepkg; Ma = 1; gam0 = 1.4; Rinner = 1; Router = 2;
%! nozzle(Ma, gam0, Rinner, Router, -1/Ma+sqrt(eps), 1, 1);

%!test
%! % Subsonic code-to-code verification against simpler nozzle1 logic
%! pkg load odepkg; Ma = 2; gam0 = 1.4; Rinner = 1.5; Router = 2;
%! [r, u, rho, p, a2, up, pp] = nozzle1(Ma, gam0, Router, Rinner, -2/7, 0.9, 1.1);
%! expected = [r(end), u(end), rho(end), p(end), a2(end), up(end), pp(end)];
%! [r, u, rho, p, a2, up, pp] = nozzle (Ma, gam0, Router, Rinner, -2/7, 0.9, 1.1);
%! observed = [r(end), u(end), rho(end), p(end), a2(end), up(end), pp(end)];
%! assert(norm(expected - observed) / norm(expected) < sqrt(sqrt(eps)))
