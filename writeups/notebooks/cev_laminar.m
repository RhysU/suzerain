#!/usr/bin/env octave

% Load a collection of problems from disk using Octave-Forge dataframe
pkg load dataframe
d = dataframe('cev_laminar.in');

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

% Save input with appended results into a new CSV-with-header file
f = fopen('cev_laminar.out','w');
dcols = columns(d);
for j = 1:dcols
  fprintf(f, '%s,', d.colnames(j))
end
scols = columns(s);
for j = 1:scols
  fprintf(f, '%s%s', s.colnames(j), merge(j != scols, ',', "\n"))
end
for i = 1:rows(d)
  t      = strrep(mat2str(d.array(i,:)), ' ', ',');
  t(end) = ',';
  fputs(f, t(2:end));
  t = strrep(mat2str(s.array(i,:)), ' ', ',');
  t(end) = "\n";
  fputs(f, t(2:end));
end
fclose(f);
