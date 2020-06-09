function r = fSimulate_response_prepare_coef_v2(nGridSpace, nFunctions)
% this function generates nFunctions functions in the spatial domain (nGridSpace)

% pre-define the smoothness
smoothness = 2;

% random function here in spatial domain
r = nan(nFunctions,nGridSpace);
for ii=1:nFunctions
%     r(ii,:) = fRandomFunction_Fourier(nGridSpace,smoothness,1);
    r(ii,:) = 10* randn(nGridSpace,1);
end

end