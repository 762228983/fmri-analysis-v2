function [logP, sla_var] = random_effects(fla_stat)

% logP = random_effects(fla_stat)
% 
% Converts first-level statistics and variance into second-level variances
% assuming random effects.
% 
% 2016-07-11: Created, Sam NH
% 
% 2016-09-21: Modified to handle NaN input, Sam NH

% second level mean
sla_mean = nanmean(fla_stat,1);

% second level variance
sla_df = sum(~isnan(fla_stat),1)-1;
sla_var = var(fla_stat,1,1) ./ sla_df;
sla_var(sla_df <= 0) = NaN;

% convert to t-stat and then p-value
logP = nan(size(sla_mean));
xi = ~isnan(sla_var);
logP(xi) = t2logP(sla_mean(xi) ./ sqrt(sla_var(xi)), sla_df(xi));

% remove first dimension (which is singleton)
dims = size(logP);
logP = reshape(logP, dims(2:end));
sla_var = reshape(sla_var, dims(2:end));