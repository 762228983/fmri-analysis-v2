function [logP, sla_var] = random_effects(fla_stat)

% logP = random_effects(fla_stat)
% 
% Converts first-level statistics and variance into second-level variances
% assuming random effects.
% 
% 2016-07-11: Created by Sam NH

% second level mean and variance
sla_mean = mean(fla_stat);
sla_df = size(fla_stat,1)-1;
sla_var = var(fla_stat,1) / sla_df;

% convert to t-stat and then p-value
logP = t2logP(sla_mean ./ sqrt(sla_var), sla_df);

% remove first dimension (which is singleton)
dims = size(logP);
logP = reshape(logP, dims(2:end));
sla_var = reshape(sla_var, dims(2:end));