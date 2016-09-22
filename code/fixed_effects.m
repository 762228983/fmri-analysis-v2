function [logP, sla_var] = fixed_effects(fla_stat, fla_var, df)

% logP = fixed_effects(fla_stat, fla_var, df)
% 
% Converts first-level statistics and variance into second-level variances
% assuming fixed effects.
% 
% 2016-07-11: Created by Sam NH
% 
% 2016-09-21: Modified to handle NaN input, Sam NH

% check that the size of fla_stat and fla_var are the same
% and that they have NaN values in the same places
assert(all(size(fla_stat)==size(fla_var)));
assert(all(isnan(fla_stat(:)) == isnan(fla_var(:))));

% expaned degrees of freedom into a matrix of size equal to fla_stat and fla_var
% and set NaN entries to have value zero
dims = size(fla_stat);
df_expanded = repmat(df, [1,dims(2:end)]);
assert(all(size(df_expanded)==size(fla_stat)));
df_expanded(isnan(fla_stat)) = 0;

% second level mean and variance
sla_mean = nanmean(fla_stat,1);
sla_var = nanmean(fla_var,1) / size(fla_stat,1);
sla_df = sum(df_expanded);

% convert to p-value
xi = sla_df > 0;
logP = nan(size(sla_mean));
logP(xi) = t2logP(sla_mean(xi) ./ sqrt(sla_var(xi)), sla_df(xi));

% remove first dimension (which is singleton)
dims = size(logP);
logP = reshape(logP, dims(2:end));
sla_var = reshape(sla_var, dims(2:end));