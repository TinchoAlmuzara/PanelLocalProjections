function LP_out = panel_LP(model)
% PANEL_LP: Panel local projection estimates and confidence intervals.
% 
% This function estimates impulse responses using local projections (LP)
% in a panel data setting and provides confidence intervals recommended in
% "Micro Responses to Macro Shocks" (Almuzara, Sancibrian). It supports
% small-sample refinements, lagged controls, and cumulative responses.
%
% USAGE:
%   LP_out = panel_LP(model)
%
% INPUT:
%   model : Structure containing the following fields:
%       - y            : (N_OBS x 1)   Regressand vector.
%       - s            : (N_OBS x N_S) Individual characteristics matrix to be interacted with shocks. Defaults to a vector of ones.
%       - X            : (N_OBS x N_X) Shock variable(s) matrix.
%       - W            : (N_OBS x N_W, optional) Control variables matrix.
%       - FE           : (N_OBS x G, optional) Fixed effects categorical matrix.
%       - i_index      : (N_OBS x 1) Individual identifiers.
%       - t_index      : (N_OBS x 1) Time period identifiers.
%       - H            : (scalar, optional) Maximum horizon. Default: ceil(0.25*T).
%       - p_max        : (scalar, optional) Maximum number of lags. Default: ceil((T-H)^(1/3)).
%       - small_sample : (logical, optional) Apply small-sample refinement. Default: false.
%       - cumulative   : (logical, optional) Compute cumulative responses. Default: false.
%
% OUTPUT:
%   LP_out : Structure containing the following fields:
%       - estimate : ((H+1) x S) Matrix of LP estimates for each horizon.
%       - SE       : ((H+1) x S) Standard errors.
%       - df       : ((H+1) x N_S) Degrees of freedom for inference (one per s-component; IK refinement gives component-specific values).
%       - CI90     : ((H+1) x S x 2) 90% confidence intervals.
%       - CI95     : ((H+1) x S x 2) 95% confidence intervals.
%       - CI99     : ((H+1) x S x 2) 99% confidence intervals.
%       - pval     : ((H+1) x S) p-values of two-sided significance tests.
%
% DEPENDENCIES:
%   - time_shift.m: Used to shift time series.
%   - regress_HDFE.m: Performs high-dimensional fixed effects regression.
%
% NOTES:
%   - t_index is either an ordered sequence of numbers or a datetime array in
%   units no smaller than a month (i.e., monthly/quarterly/annual/biennial
%   units are allowed); t_index is rescaled so that the minimum time-step
%   is one unit, and observations for which t_index is not a multiple of
%   the minimum time-step are treated as missing.
%   - The horizon H and the lag length p_max are measured in the normalized
%   unit of t_index.
%   - Small-sample refinements follow Imbens & Kolesar (2016, REStat), with
%   asymptotic time-clustered standard errors used if small_sample = false.
%
% Version: 2024 Jun 10

% Extract inputs
y       = model.y;
if isfield(model, 's'); s = model.s; else, s = ones(size(y, 1), 1); end
X       = model.X;
i_index = model.i_index;
t_index = model.t_index;

% Normalize time indexes if needed
if (isdatetime(t_index) == true)
    t_index = 12*(year(t_index) - min(year(t_index))) + month(t_index);
end
t_min   = min(t_index);
t_diff  = min(diff(sort(unique(t_index), 'ascend')));
t_index = (t_index - t_min)/t_diff + 1;

% Drop units with inconsistent time indexes
keep    = ~(t_index - floor(t_index) > 1e-2);
y       = y(keep);
s       = s(keep, :);
X       = X(keep, :);
i_index = i_index(keep);
t_index = round(t_index(keep));

% Recover dimensions
n_obs = length(t_index);
T_eff = length(unique(t_index));
n_s   = size(s, 2)*size(X, 2); 

% Construct interacted regressor
sX = reshape( repmat(s, [1, 1, size(X, 2)]) .* permute(repmat(X, [1, 1, size(s, 2)]), [1, 3, 2]), [n_obs, n_s]);

% Recover optional data
if isfield(model, 'W')  && ~isempty(model.W),  W  = model.W(keep, :);  else, W  = zeros(n_obs, 0); end
if isfield(model, 'FE') && ~isempty(model.FE), FE = model.FE(keep, :); else, FE = zeros(n_obs, 0); end

% Recover optional arguments
if isfield(model, 'H'),            H            = model.H;            else, H            = ceil(0.25*T_eff);      end
if isfield(model, 'p_max'),        p_max        = model.p_max;        else, p_max        = ceil((T_eff-H)^(1/3)); end
if isfield(model, 'small_sample'), small_sample = model.small_sample; else, small_sample = false;                 end
if isfield(model, 'cumulative'),   cumulative   = model.cumulative;   else, cumulative   = false;                 end

% Preallocate output
LP_estimate = zeros(H+1, n_s);
LP_SE       = zeros(H+1, n_s);
LP_df       = zeros(H+1, n_s);
LP_CI90     = zeros(H+1, n_s, 2);
LP_CI95     = zeros(H+1, n_s, 2);
LP_CI99     = zeros(H+1, n_s, 2);
LP_pval     = zeros(H+1, n_s);

% Iterate over horizons
y_h = zeros(n_obs, 1);
for h = 0:H

    % Lead regressand
    if (cumulative == true)
        y_h = y_h + time_shift(y, i_index, t_index, h);
    else
        y_h = time_shift(y, i_index, t_index, h);
    end

    % Construct lagged controls
    p     = min(h, p_max);
    W_lag = NaN(n_obs, p, 1+n_s);
    for j = 1:p
        W_lag(:, j, 1) = time_shift(y, i_index, t_index, -j);
        for i_s = 1:n_s
            W_lag(:, j, 1+i_s) = time_shift(sX(:, i_s), i_index, t_index, -j);
        end
    end
    W_lag = reshape(W_lag, [length(i_index), (1+n_s)*p]);

    % Prepare data
    d      = ~any(isnan([y_h, sX, W_lag, W, FE]), 2);
    y_LP   = y_h(d, :);
    X_LP   = [sX(d, :), W_lag(d, :), W(d, :)];
    n_X    = size(X_LP, 2);
    dum_LP = FE(d, :);

    % Prepare time-series indexes
    t_LP  = t_index(d);
    t_set = sort(unique(t_LP), 'ascend');
    T     = length(t_set);
    
    % Compute LP estimator    
    [b_LP, y_LP, X_LP]  = regress_HDFE(y_LP, X_LP, dum_LP);
    LP_estimate(h+1, :) = b_LP(1:n_s);

    % Compute score and hessian
    Xv_it = X_LP .* (y_LP - X_LP*b_LP);
    Xv_t  = zeros(T, n_X);
    for t = 1:T
        t_tmp      = (t_LP == t_set(t));
        Xv_t(t, :) = sum(Xv_it(t_tmp, :), 1);
    end
    XX    = (X_LP')*X_LP;

    % Compute t-LAHR standard error
    if (small_sample == true)

        for i_s = 1:n_s
            % Compute Imbens-Kolesar small-sample refinement
            s_LP = s(d, i_s);
            X_t  = zeros(T, n_X);
            ss_t = zeros(T, 1);                       % per-period s'*s denominator
            for t = 1:T
                t_tmp     = (t_LP == t_set(t));
                ss_t(t)   = s_LP(t_tmp)'*s_LP(t_tmp);
                X_t(t, :) = (s_LP(t_tmp)'*X_LP(t_tmp, :)) / ss_t(t); % ill-conditioning checked below
            end
            warn_if_illconditioned_ik(ss_t, X_t, h, i_s);
            % A truly-degenerate period (s'*s = 0) makes X_t non-finite. Set df/SE to
            % NaN for this cell and skip so a single degenerate period does not abort the
            % run (rather than being silently reported as SE = 0).
            if (~all(isfinite(X_t(:))))
                LP_SE(h+1, i_s) = NaN;
                LP_df(h+1, i_s) = NaN;
                continue
            end
            P0     = eye(T) - X_t*pinv((X_t')*X_t)*(X_t');
            Xv_var = ((Xv_t./sqrt(diag(P0)))')*(Xv_t./sqrt(diag(P0)));
            b_var  = pinv(XX)*Xv_var*pinv(XX);
            G0     = zeros(T);
            XX0    = pinv((X_t')*X_t);
            for t = 1:T, G0(:, t) = P0(:, t)*X_t(t, :)*XX0(:, 1)/sqrt(P0(t, t)); end
            lam0   = eig(G0'*G0);

            % Store standard error and degrees of freedom
            LP_SE(h+1, i_s) = sqrt(max(0, diag(b_var(i_s, i_s))));
            LP_df(h+1, i_s) = (sum(lam0))^2/(sum(lam0.^2));
        end

    else

        % Compute time-clustered sandwich formula
        Xv_var = (Xv_t')*Xv_t;
        b_var  = pinv(XX)*Xv_var*pinv(XX);

        % Store standard error and degrees of freedom
        LP_SE(h+1, :) = sqrt(max(0, diag(b_var(1:n_s, 1:n_s))));
        LP_df(h+1, :) = Inf;

    end

    % Compute confidence intervals (LP_df(h+1, :) is a 1×n_s row vector;
    % tinv broadcasts element-wise so each s-component uses its own df)
    cv_tmp             = tinv(1-(1-0.90)/2, LP_df(h+1, :));
    LP_CI90(h+1, :, 1) = LP_estimate(h+1, :) - cv_tmp .* LP_SE(h+1, :);
    LP_CI90(h+1, :, 2) = LP_estimate(h+1, :) + cv_tmp .* LP_SE(h+1, :);
    cv_tmp             = tinv(1-(1-0.95)/2, LP_df(h+1, :));
    LP_CI95(h+1, :, 1) = LP_estimate(h+1, :) - cv_tmp .* LP_SE(h+1, :);
    LP_CI95(h+1, :, 2) = LP_estimate(h+1, :) + cv_tmp .* LP_SE(h+1, :);
    cv_tmp             = tinv(1-(1-0.99)/2, LP_df(h+1, :));
    LP_CI99(h+1, :, 1) = LP_estimate(h+1, :) - cv_tmp .* LP_SE(h+1, :);
    LP_CI99(h+1, :, 2) = LP_estimate(h+1, :) + cv_tmp .* LP_SE(h+1, :);

    % Compute p-values of two-sided significance test
    LP_pval(h+1, :) = 2*(1-tcdf(abs(LP_estimate(h+1, :) ./ LP_SE(h+1, :)), LP_df(h+1, :)));

end

% Store output
LP_out          = struct();
LP_out.estimate = LP_estimate;
LP_out.SE       = LP_SE;
LP_out.df       = LP_df;
LP_out.CI90     = LP_CI90;
LP_out.CI95     = LP_CI95;
LP_out.CI99     = LP_CI99;
LP_out.pval     = LP_pval;

end


%% LOCAL FUNCTIONS

function y1 = time_shift(y0, i_index, t_index, L)
% TIME_SHIFT: Compute lag or lead of time series.
%
% Inputs:
%   Y0      - An n-by-m matrix of time series data, where n is the number of observations
%             and m is the number of variables.
%   i_index - An n-by-1 vector indicating the unit (individual) index for each observation.
%   t_index - An n-by-1 vector indicating the time index for each observation.
%   L       - An integer indicating the lead (if positive) or lag (if negative)
%             to be applied to the time series data.
%
% Outputs:
%   Y1 - An n-by-m matrix of time-shifted data, where each element is the lead or lag
%        of the corresponding element in y0, as specified by L.
%
% Version: 2024 Jun 10

% Loop over units
y1 = NaN(size(y0));
for i = unique(i_index)'

    % Extract unit data
    i_unit  = (i_index == i);
    t_unit  = t_index(i_unit);
    y0_unit = y0(i_unit, :);

    % Obtain data lag (negative L) or lead (positive L)
    [is_t, loc_t]    = ismember(t_unit+L, t_unit);
    y1_unit          = NaN(size(y0_unit));
    y1_unit(is_t, :) = y0_unit(loc_t(is_t), :);
    y1(i_unit, :)    = y1_unit;

end

end

function [b, y_resid, X_resid] = regress_HDFE(y, X, FE, tol, max_iter)
% REGRESS_HDFE: Efficient high-dimensional fixed effects regression.
%
% INPUTS:
%   y        - (N x 1) Dependent variable
%   X        - (N x K) Independent variables
%   FE       - (N x G) Group identifiers (e.g., {firm, time, industry})
%   tol      - Convergence tolerance (default: 1e-8)
%   max_iter - Maximum number of iterations (default: 500)
%
% OUTPUTS:
%   b       - (K x 1) Estimated coefficients
%   Y_RESID - (N x 1) Dependent variable orthogonalized wrt fixed effects
%   X_RESID - (N x K) Independent variables orthogonalized wrt fixed effects
%
% Version: 2024 Jun 10

% Set default options
if nargin < 4, tol      = 1e-8; end
if nargin < 5, max_iter = 500;  end

% Recover number of regressors
n_X = size(X, 2);

% Convert identifiers to unique indices
n_FE      = size(FE, 2);
group_ids = cell(1, n_FE);
for i_FE = 1:n_FE
    [~, ~, group_ids{i_FE}] = unique(FE(:, i_FE));
end

% Initialize transformed variables
y_resid   = y;
X_resid   = X;
converged = false;
iter      = 0;

% Perform iterative residualization
while ~converged && (iter < max_iter)
    iter  = iter + 1;
    X_old = X_resid;
    y_old = y_resid;
    
    for i_FE = 1:n_FE
        % Compute means by fixed effect group
        group_mean_y = accumarray(group_ids{i_FE}, y_resid, [], @mean);
        group_mean_X = zeros(length(group_mean_y), n_X);
        for i_X = 1:n_X
            group_mean_X(:, i_X) = accumarray(group_ids{i_FE}, X_resid(:, i_X), [], @mean);
        end

        % Remove fixed effect from y and X
        y_resid = y_resid - group_mean_y(group_ids{i_FE});
        if (n_X > 0)
            X_resid = X_resid - group_mean_X(group_ids{i_FE}, :);
        end
    end
    
    % Check convergence
    if norm([y_resid, X_resid] - [y_old, X_old], 'fro') < tol
        converged = true;
    end
end

if ~converged
    warning('Fixed effects did not converge after %d iterations.', max_iter);
end

% Singularity safeguard: emit a STRONG warning when the residualized design matrix
% is severely ill-conditioned -- i.e. a regressor is (near-)collinear with the
% others, so the specification is effectively unidentified. At high condition
% numbers the OLS coefficient estimates are numerically unreliable. cond() is the
% 2-norm condition number. Detect-and-report only: it does not alter the numerics
% and is a no-op on well-conditioned specs.
if (n_X > 0)
    kappa_X = cond(X_resid);
    if kappa_X > 1e5
        warning(['SINGULARITY SAFEGUARD: the design matrix is severely ill-conditioned ' ...
                 '(condition number = %.3g > 1e5). One or more regressors are nearly ' ...
                 'collinear with the rest, so the specification is effectively ' ...
                 'unidentified: the OLS coefficient estimates are numerically unreliable. ' ...
                 'Check the specification for a redundant or near-duplicate control / ' ...
                 'interaction term.'], kappa_X);
    end
end

% Run OLS on residualized data
if (n_X > 0)
    b = X_resid\y_resid;
else
    b = [];
end

end

% Helper: warn_if_illconditioned_ik
% IK small-sample (X_t) ill-conditioning detector. WARNs that a specification is
% ill-conditioned rather than silently rounding it off with a numerical floor.
% Fires on EITHER trigger, per (h, s-component):
%   (1) a within-period instrument sum-of-squares s'*s below ss_tol (where the
%       per-period projection denominator is ~0 -- the projection is ~0/0); or
%   (2) the 2-norm condition number of the IK design X_t exceeding kappa_max.
% DETECT-ONLY; never alters the numerics.
function warn_if_illconditioned_ik(ss_t, X_t, h, i_s, ss_tol, kappa_max)
    if (nargin < 5 || isempty(ss_tol)),    ss_tol    = 1e-12; end
    if (nargin < 6 || isempty(kappa_max)), kappa_max = 1e8;   end

    ss_min = min(ss_t);
    if (~isfinite(ss_min) || ss_min < ss_tol)
        warning('panelLP:IKillConditioned', ...
            ['ILL-CONDITIONED IK SPEC (h=%d, s-component %d): a within-period instrument ' ...
             'sum-of-squares s''s = %.3g fell below tolerance %.0g. The per-period projection ' ...
             'denominator is ~0, so the small-sample (X_t) refinement is unstable and the df/SE ' ...
             'for this cell are unreliable. Consider rescaling the instrument or dropping the ' ...
             'degenerate period.'], ...
            h, i_s, ss_min, ss_tol);
    end

    if (all(isfinite(X_t(:))))
        kappa_X = cond(X_t);
    else
        kappa_X = Inf;
    end
    if (kappa_X > kappa_max)
        warning('panelLP:IKillConditioned', ...
            ['ILL-CONDITIONED IK SPEC (h=%d, s-component %d): the small-sample design matrix X_t ' ...
             'is severely ill-conditioned (condition number = %.3g > %.0g). The Imbens-Kolesar ' ...
             'degrees-of-freedom estimate is numerically unstable; treat the small-sample ' ...
             'inference for this cell with caution.'], ...
            h, i_s, kappa_X, kappa_max);
    end
end
