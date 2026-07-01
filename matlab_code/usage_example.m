%%% USAGE EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script illustrates the usage of panel_LP, the function to implement
% estimation and inference as recommended in "Micro Responses to Macro Shocks" 
% by M. Almuzara and V. Sancibrian.
%
% Version: 2025 February 19 - Matlab R2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear memory
clear
close all
clc
rng(1916)

% Simulate a simple dataset
T       = 30;
N       = 1000;
X       = repmat(randn(1, T), [N, 1]);
X       = X(:);
s       = repmat(1 + randn(N, 1), [1, T]);
s       = s(:);
b_true  = 1;
y       = (s .* X) * b_true + reshape((0.5*s(1:N) + rand(N, 1)).*randn(1, T), [N*T, 1]) + randn(N*T, 1);
t_index = repmat((1:T), [N, 1]);
t_index = t_index(:);
i_index = repmat((1:N)', [1, T]);
i_index = i_index(:);

% Create model structure to pass inputs
model              = struct();
model.y            = y;
model.s            = s;
model.X            = X;
model.i_index      = i_index;
model.t_index      = t_index;
model.W            = [];                 % this may be omitted if no controls are to be used
model.FE           = [i_index, t_index]; % unit and time fixed effects, may be omitted if no fixed effects are to be used
model.H            = 5;                  % impulse response horizon
model.p_max        = 0;                  % number of lags of regressand and shock to be added as controls
model.small_sample = true;               % Imbens-Kolesar-2016-REStat small sample refinement; defaults to false
model.cumulative   = true;               % Report cumulative impulse responses; default to false

% Call panel local projections function
LP_out = panel_LP(model);

% Plot estimates with 90% confidence bands
H        = model.H;
horizons = (0:H)';
fig0  = figure();
ax0   = axes();
hold on
fill(ax0, [horizons; flipud(horizons)], [LP_out.CI90(:, 1, 1); flipud(LP_out.CI90(:, 1, 2))], [  0,  85, 164]/255, 'facealpha', 0.2, 'edgecolor', 'none')
plot0 = plot(ax0, horizons, [LP_out.CI90(:, :, 1), LP_out.estimate, LP_out.CI90(:, :, 2)]);
hold off
xlabel(ax0, 'h')
ylabel(ax0, 'LP estimates of cumulative impulse responses');

% Tune ax handle
set(ax0, 'LooseInset', get(ax0, 'TightInset'))
xlim(ax0, [0, H])
xticks(ax0, horizons)
ylim(ax0, [-0.7, 2.7])
yticks(ax0, -0.5:0.5:2.5)
set(ax0, 'fontname', 'helvetica', 'fontsize', 25)
grid(ax0, 'on')
line(ax0, [0, H], [0, 0], 'linewidth', 1, 'color', [0, 0, 0], 'linestyle', '--')

% Tune plot
legend(ax0, {'90% CI', '', 'Point estimate', ''}, 'location', 'best', 'box', 'off')
set(plot0, {'linewidth'}, {4; 6; 4})
set(plot0, 'color', [  0,  85, 164]/255)
set(plot0, {'linestyle'}, {':'; '-'; ':'})
set(plot0, {'marker'}, {'none'; 'o'; 'none'})
set(plot0, 'markersize', 15)

% Tune figure
set(fig0, 'theme', 'light', 'units', 'inches', 'position', [0 0 14 14], 'paperunits', 'inches', 'papersize', [14 14], 'paperposition', [0 0 14 14])
