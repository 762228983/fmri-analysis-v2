function X_hrf_convolved =  convolve_with_hrf(X, hrf_name, sr)

% Convolves a design matrix X, with an hrf

% get hrf
switch hrf_name
    case 'fsfast-gamma-BOLD'
        hrf = hrf_fsfast_gamma(1/sr, 'BOLD', 'noplot');
        hrf = hrf/abs(sum(hrf));
    otherwise
        error('hrf_name %s not valid\n', hrf_name);
end

% pad (typical) or truncate (unusual) hrf
N = size(X,1);
if length(hrf) < N
    hrf = [hrf; zeros(N - length(hrf),1)];
elseif length(hrf) > N
    hrf = hrf(1:N);
end

% convolve design matrix with hrf
FT_X = fft(X);
FT_hrf = fft(hrf);
X_hrf_convolved = ...
    ifft(FT_X .* repmat(FT_hrf, 1, size(FT_X,2)));