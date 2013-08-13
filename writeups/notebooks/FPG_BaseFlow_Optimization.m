% Script file loading necessary packages for use within Octave
1; if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end

function h = baseflow_h(x)
% Compute nonlinear constraints for use by baseflow_sqp.  See baseflow_sqp.
% Specify h(x) >= 0 constraints using eps to accomplish h(x) > 0 where needed.
% For subsonic, -1/Ma < u1 < 0 achieves at most nearly sonic outflow.
% For supersonic, u1 > 1/Ma achieves at least nearly sonic inflow.
% Always R0 >= R1 must hold to maintain the orientation of the domain.
  [dp_e, dstar, gam0, Ma, Ma_e, p1, R0, R1, rho1, T_e, u1] = num2cell(x){:};
  if Ma_e < 1
      h = [R0-R1; 1/Ma-u1-eps; u1-eps];
  else
      h = [R0-R1; u1-1/Ma-eps];
  end
end

function phi = baseflow_phi(x)
% Compute nonlinear functional for use by baseflow_sqp.  See baseflow_sqp.
% Compute phi(x) returning the l_1 of the absolute mismatch in dp_e, Ma_e, T_e.
% The radial problem is solved by nozzle(...) on the smallest possible domain.
  [dp_e, dstar, gam0, Ma, Ma_e, p1, R0, R1, rho1, T_e, u1] = num2cell(x){:};
  R2 = sqrt(R0**2 + dstar**2);
  [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1);
  phi = abs(Ma_e - Ma*r(end)*abs(u(end)) / (R2*sqrt(a2(end)))) ...
      + abs(dp_e + R2*abs(pp(end)) / R0                      ) ...
      + abs(T_e  - a2(end)                                   );
end

function s = baseflow_sqp(dp_e, dstar, gam0, Ma_e, T_e)
% Driver for achieving target dp_e, Ma_e, and T_e quantities at (R0, dstar).
% Establish a sequential quadratic programming problem and solve with sqp:
%    minimize phi(x) subject to g(x) = 0, h(x) >= 0, l <= x <= u
% where x is the (alphabetically-ordered) vector which may be unpacked with
%    [dp_e, dstar, gam0, Ma, Ma_e, p1, R0, R1, rho1, T_e, u1] = num2cell(x){:};
% and l, u are used to hold some parameters fixed.  See documentation for
% http://www.gnu.org/software/octave/doc/interpreter/Nonlinear-Programming.html.

  % Establish initial guess, lower bounds, and upper bounds for each parameter
  x    =zeros(11, 1); l    =zeros(size(x)); u    =zeros(size(x));
  x( 1)=dp_e;         l( 1)=dp_e;           u( 1)=dp_e;     % dp_e  constant
  x( 2)=dstar;        l( 2)=dstar;          u( 2)=dstar;    % dstar constant
  x( 3)=gam0;         l( 3)=gam0;           u( 3)=gam0;     % gam0  constant
  x( 4)=1;            l( 4)=eps;            u( 4)=realmax;  % Ma    positive
  x( 5)=Ma_e;         l( 5)=Ma_e;           u( 5)=Ma_e;     % Ma_e  constant
  x( 6)=1;            l( 6)=eps;            u( 6)=realmax;  % p1    positive
  x( 7)=1;            l( 7)=eps;            u( 7)=realmax;  % R0    positive
  x( 8)=1;            l( 8)=eps;            u( 8)=realmax;  % R1    positive
  x( 9)=1;            l( 9)=eps;            u( 9)=realmax;  % rho1  positive
  x(10)=T_e;          l(10)=T_e;            u(10)=T_e;      % T_e   constant
  x(11)=NaN;          l(11)=-realmax;       u(11)=realmax;  % u1    unleashd

  % Guess for u1 = x(11) must be consistent with sub- vs supersonic conditions
  x(11) = sign(Ma_e - 1) * Ma_e + eps;

  % Run the sequential quadratic programming algorithm provided by Octave
  % placing results into a struct.  On success, s.info == 101 and s.nozzle is
  % curried so that s.nozzle(Ly) provides data on (R0,0) to (R0,Ly).
  s = struct('x0', x);
  [s.x, s.obj, s.info, s.iter, s.nf, s.lambda]                         ...
        = sqp(x, @baseflow_phi, [], @baseflow_h, l, u, 1);
  [s.dp_e,s.dstar,s.gam0,s.Ma,s.Ma_e,s.p1,s.R0,s.R1,s.rho1,s.T_e,s.u1] ...
        = num2cell(s.x){:};
  s.nozzle=@(Ly) nozzle(s.Ma,s.gam0,s.R1,sqrt(s.R0**2+Ly**2),s.u1,s.rho1,s.p1);

end

% A sample invocation in the spirit of 4m leeward of the stagnation point
% s = baseflow_sqp(-8.75, 1, 1.4, 1.1, 4.2)
