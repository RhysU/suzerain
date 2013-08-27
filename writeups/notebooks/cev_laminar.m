#!/usr/bin/env octave

% Load a collection of problems from disk using Octave-Forge dataframe
pkg load dataframe
d = dataframe('cev_laminar.csv');

% Process problems in parallel solving with delta fixed to be one
pkg load general
tic()
s = parcellfun(max(1, nproc()-1), @nozzle_baseflow,      ...
               num2cell(ones(rows(d),1)),                ...
               num2cell(d.gamma_e), num2cell(d.Ma_edge), ...
               num2cell(d.p_exi), num2cell(d.T_ratio));
toc()

% Convert the returned results into a dataframe
s = dataframe(struct2cell(s).', 'colnames', fieldnames(s));
