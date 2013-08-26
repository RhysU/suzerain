#!/usr/bin/env octave
d = dlmread('cev_laminar.csv', ',', 1, 0);
clear s;
for i = 1:rows(d)
    tic
    try
        s(i) = nozzle_baseflow(1.0, d(i,3), d(i,4), d(i,5), d(i,6),
                               optimset('MaxIter', 100,
                                        'TolFun', sqrt(eps),
                                        'weights', [1;1;1]));
        s(end)
    catch
        s(i).cvg = 0
        warning('nozzle_baseflow(%g, %g, %g, %g, %g) fails', ...
                d(i,2), d(i,3), d(i,4), d(i,5), d(i,6));
    end
    toc
end
