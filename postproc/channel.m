#!/usr/bin/octave -qf
% A GNU Octave script for exploratory post-processing on the HDF5 output of
% Suzerain's channel_mean.  For example, after aggregating samples
%
%    shell$    channel_mean -o samples.h5 sample*.h5
%
% open GNU Octave, load the collection of samples, and run this script
% to batch process the statistics into an explorable form:
%
%    octave:1> load samples.h5
%    octave:2> channel
%
% Running 'whos' will give an indication of the available statistics.  They
% statistics are named and described per the output of 'channel_mean -d'.
% For example, 'bar_u' is the Reynolds-averaged velocity and 'tilde_nu'
% is the Favre-averaged kinematic viscosity.  Running this script twice
% will result in errors.  If you have plots or processing you'd like to
% repeatedly perform, place it in its own script to be run afterwards
%
%    octave:3> after_channel_the_flood
%
% Each statistic is stored in an Octave struct with fields named according
% to 'help arsel'.  For example, 'bar_u.mu' is the mean velocity and
% 'bar_u.mu_sigma' its estimated standard deviation.  Two dimensional matrices
% of raw samples on the space-time grid 'y' by 't' are available as, e.g.,
% 'bar_u.data'.  Generally you want to compute derived quantities as pointwise
% information from '*.data' and to then run this data through the function
% arsel(...).  See 'tau_u' processing below as an example.
format compact;
more off;

% Ensure arsel/arcov post-processing utilities are available.
% Build these somewhere and use addpath(...) to inform Octave they exist.
if exist('arsel', 'file') != 3
    error('Did not detect arsel(...).  See https://github.com/rhysu/ar');
end
if exist('arcov', 'file') != 3
    error('Did not detect arcov(...).  See https://github.com/rhysu/ar');
end

% Load specified files arguments if invoked from command line
% TODO Permit batch processing across multiple data sets
if strcmp(program_invocation_name(),"octave")
    prefix = '';
else
    for i = 1:nargin
        printf('Loading %s\n', argv(){i});
        load(argv(){i});
        [dir, prefix, ext] = fileparts(argv(){i});
        clear dir ext
    end
end

printf('Processing all available samples in the workspace\n');
tic();
S = cat(1, who('bar_*'),   ...            % Reynolds averages
           who('tilde_*'), ...            % Favre averages
           who('local_*'));               % Local quantities
S = strcat(S, {' = arsel('}, S, {');'});  % Build in-place arsel commands
parfor i = 1:length(S)
    printf('Processing %d/%d: %s\n', i, length(S), cell2mat(S(i)));
    eval(cell2mat(S(i)));
end
printf('Processing completed in %g seconds\n', toc());

tic();
printf('Computed classical nondimensional quantities');
% Notice scaling factors required due to nondimensionalization
tau_w    = arsel(bar_mu.data(1,:).*bar_u__y.data(1,:));
u_tau    = arsel(sqrt(1/Re)*sqrt(tau_w.data ./ bar_rho.data(1,:)));
delta_nu = arsel((1/Re) * tilde_nu.data(1,:) ./ u_tau.data);
Re_tau   = arsel((max(y) - min(y)) ./ (2 .* delta_nu.data));
eta      = arsel((Re^(-3) * tilde_nu.data.**3 ./ tilde_epsilon.data).**(0.25));
printf(' in %g seconds\n', toc());

printf('Displaying basic grid metrics\n');
Nx
Ny
Nz
Lx
Ly
Lz
delta_x = Lx / double(Nx)
delta_z = Lz / double(Nz)
% Everything in plus units "lives" within the plus struct
plus   = struct('delta_x', delta_x / delta_nu.mu, ...
                'delta_z', delta_z / delta_nu.mu)
plus.y = y / delta_nu.mu;
printf('Collocation points below plus.y of 10:\n');
disp(plus.y(plus.y <= 10));

tic();
printf('Computed TKE-related quantities in plus units');
plus.upp_upp = arsel(tilde_upp_upp.mu ./ u_tau.mu**2);
plus.vpp_vpp = arsel(tilde_vpp_vpp.mu ./ u_tau.mu**2);
plus.wpp_wpp = arsel(tilde_wpp_wpp.mu ./ u_tau.mu**2);
plus.upp_vpp = arsel(tilde_upp_vpp.mu ./ u_tau.mu**2);
plus.k       = arsel(tilde_k      .mu ./ u_tau.mu**2);
printf(' in %g seconds\n', toc());

%% TODO KRUFT that should be moved elsewhere
%% Plot TKE-related sample means
%h = figure;
%set(h,'papertype', '<custom>')
%set(h,'paperunits','inches');
%set(h,'papersize',[3.5 7])
%set(h,'paperposition', [0,0,[7 3.5]])
%subplot(1,2,1);
%plot(plus.y, plus.upp_upp.mu, ...
%     plus.y, plus.vpp_vpp.mu, ...
%     plus.y, plus.wpp_wpp.mu, ...
%     plus.y, plus.upp_vpp.mu, ...
%     plus.y, plus.k      .mu);
%legend('u''''u''''', ...
%       'v''''v''''', ...
%       'w''''w''''', ...
%       'u''''v''''', ...
%       'k');
%xlim([0, round(Re_tau.mu)]);
%ylabel('Normalized by utau2');
%xlabel('yplus');
%subplot(1,2,2);
%plot(plus.y, plus.upp_upp.mu_sigma, ...
%     plus.y, plus.vpp_vpp.mu_sigma, ...
%     plus.y, plus.wpp_wpp.mu_sigma, ...
%     plus.y, plus.upp_vpp.mu_sigma, ...
%     plus.y, plus.k      .mu_sigma);
%legend('std u''''u''''', ...
%       'std v''''v''''', ...
%       'std w''''w''''', ...
%       'std u''''v''''', ...
%       'std k');
%xlim([0, round(Re_tau.mu)]);
%xlabel('yplus');
%ylabel('Stddev normalized by utau2');
%set(findall(h,'-property','LineWidth'), 'LineWidth', 4);
%subplot(1,2,1);
%set(gca, 'LineWidth', 2);
%set(gca, 'Position', [0.10000 0.11000 0.33466 0.81500]);
%subplot(1,2,2);
%set(gca, 'LineWidth', 2);
%set(gca, 'Position', [0.59034 0.11000 0.33466 0.81500]);
%print -depsc2 sample_means_fm.eps;
%
%% Plot u''u'' covariance structures
%M = arcov(plus.upp_upp);
%
%h = figure;
%contourf(plus.y(1:end/2), plus.y(1:end/2), M.eff_cov(1:end/2,1:end/2), 15)
%colorbar('EastOutside');
%ylabel('yplus');
%xlabel('yplus');
%title("Autocorrelation of u''u''")
%axis('square');
%print -depsc2 sample_covuu_fm.eps;
%
%h = figure;
%contourf(plus.y(1:end/2), plus.y(1:end/2), M.eff_N(1:end/2,1:end/2), 15)
%colorbar('EastOutside');
%ylabel('yplus');
%xlabel('yplus');
%title("Effectively independent samples of u''u''")
%axis('square');
%print -depsc2 sample_covuuN_fm.eps;
%
%% Plot u''v'' covariance structures
%M = arcov(plus.upp_vpp);
%
%h = figure;
%contourf(plus.y(1:end/2), plus.y(1:end/2), M.eff_cov(1:end/2,1:end/2), 15)
%colorbar('EastOutside');
%ylabel('yplus');
%xlabel('yplus');
%title("Autocorrelation of u''v''")
%axis('square');
%print -depsc2 sample_covuv_fm.eps;
%
%h = figure;
%contourf(plus.y(1:end/2), plus.y(1:end/2), M.eff_N(1:end/2,1:end/2), 15)
%colorbar('EastOutside');
%ylabel('yplus');
%xlabel('yplus');
%title("Effectively independent samples of u''v''")
%axis('square');
%print -depsc2 sample_covuvN_fm.eps;
