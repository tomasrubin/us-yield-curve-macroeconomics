function covLags = fEstimateCovLags(onb,specDensity,numOfLags_all)

covLags = [];

% integrate spectral density to get covariances
covLags.ONB = zeros(numOfLags_all+1,onb.nBasis,onb.nBasis);
for lag=0:numOfLags_all
    m = zeros(onb.nBasis);
    for k=1:specDensity.nGridFreq % Fourier the positified spectral density
        omega = specDensity.gridFreq(k);
        m = m + squeeze(specDensity.ONB_positified(k,:,:)) * exp(1i*omega*lag) / specDensity.nGridFreq*2*pi;
    end    
    covLags.ONB(lag+1,:,:) = real(m); % covariances are real, delete the numerical imag part
end

% crate shortcut for cov lag 0
covLags.lag0 = squeeze( covLags.ONB(1,:,:) );

% calculate the norms
covLags.norms = zeros(numOfLags_all+1,1);
for lag=0:numOfLags_all
    covLags.norms(lag+1) = norm( squeeze(covLags.ONB(lag+1,:,:)) );
end

% calculate the (recommended) cut-of-lag
covLags.cutOfLag = find( covLags.norms < 0.005 , 1);
if (~isnumeric(covLags.cutOfLag) || isempty(covLags.cutOfLag)), covLags.cutOfLag = 100; end
if (covLags.cutOfLag > 101), covLags.cutOfLag = 100; end
if (covLags.cutOfLag < 10), covLags.cutOfLag = 10; end


end