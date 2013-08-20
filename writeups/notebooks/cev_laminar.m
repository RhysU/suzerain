#!/usr/bin/env octave
guess = [];
d = dlmread('cev_laminar.csv', ',', 1, 0);
clear s;
for i = 1:rows(d)
    tic
    try
        s(i) = nozzle_baseflow(d(i,2), d(i,3), d(i,4), d(i,5), d(i,6),
                             eps, 100, 0, guess);
        guess = [s(i).Ma; s(i).R0; s(i).rho1; s(i).u1];
        s(end)
    catch
        warning('nozzle_baseflow(%g, %g, %g, %g, %g) fails', ...
                d(i,2), d(i,3), d(i,4), d(i,5), d(i,6));
    end
    toc
end
