function [logP, sla_var] = fixed_effects(fla_stat, fla_var, df)

% logP = fixed_effects(fla_stat, fla_var, df)
% 
% Converts first-level statistics and variance into second-level variances
% assuming fixed effects.
% 
% 2016-07-11: Created by Sam NH

% second level mean and variance
sla_mean = mean(fla_stat);
sla_var = mean(fla_var) / size(fla_stat,1);
sla_df = sum(df);

% convert to p-value
logP = t2logP(sla_mean ./ sqrt(sla_var), sla_df);

% remove first dimension (which is singleton)
dims = size(logP);
logP = reshape(logP, dims(2:end));
sla_var = reshape(sla_var, dims(2:end));


