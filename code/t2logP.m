function logP = t2logP(t, df)

assert(all(df(:) > 0));
t(isinf(t) | t == 0) = NaN;
logP = -sign(t).*log10(2*tpvalue_copy(-abs(t), df)); % two-tailed