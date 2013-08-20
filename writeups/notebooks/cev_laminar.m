#!/usr/bin/env octave
d = dlmread('cev_laminar.csv', ',', 1, 0);
guess = [];
clear s;
for i = 10:rows(d)
    tic
    s(i) = nozzle_baseflow(d(i,2), d(i,3), d(i,4), d(i,5), d(i,6),
                           sqrt(eps), 100, 1, guess);
    toc
    guess = [s(i).Ma; s(i).R0; s(i).rho1; s(i).u1];
end
