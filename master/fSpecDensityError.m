function specDensity_error = fSpecDensityError( onb, specDensity, specDensity_true )

specDensity_error = [];

%% frobenius norm
specDensity_error.mse_fro = 0;
specDensity_error.r2mse_fro = 0;
rmse_divisor_fro = 0;

for k = 1:specDensity.nGridFreq
    m = onb.onbMatrix*(squeeze(specDensity.ONB(k,:,:)) - squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    specDensity_error.mse_fro = specDensity_error.mse_fro + 2*pi/specDensity.nGridFreq * norm(m,'fro')^2 / onb.nGridSpace^2;

    m2 = onb.onbMatrix*(squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    rmse_divisor_fro = rmse_divisor_fro + 2*pi/specDensity.nGridFreq * norm(m2,'fro')^2 / onb.nGridSpace^2;
    specDensity_error.r2mse_fro = specDensity_error.r2mse_fro +...
        2*pi/specDensity.nGridFreq * (norm(m,'fro')^2) / (  norm(m2,'fro')^2 );
end

specDensity_error.rmse_fro = specDensity_error.mse_fro / rmse_divisor_fro;

%% integrate maximum norm
specDensity_error.mse_intmax = 0;
specDensity_error.r2mse_intmax = 0;
rmse_divisor_intmax = 0;

for k = 1:specDensity.nGridFreq
    m = onb.onbMatrix*(squeeze(specDensity.ONB(k,:,:)) - squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    specDensity_error.mse_intmax = specDensity_error.mse_intmax + 2*pi/specDensity.nGridFreq * max(abs(m(:)));

    m2 = onb.onbMatrix*(squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    rmse_divisor_intmax = rmse_divisor_intmax + 2*pi/specDensity.nGridFreq * max(abs(m2(:)));
    specDensity_error.r2mse_intmax = specDensity_error.r2mse_intmax +...
        2*pi/specDensity.nGridFreq * ( max(abs(m(:))) / max(abs(m2(:))) );

end

specDensity_error.rmse_intmax = specDensity_error.mse_intmax / rmse_divisor_intmax;

%% max-max norm          
m = zeros(specDensity.nGridFreq,1);
m2 = zeros(specDensity.nGridFreq,1);

for k = 1:specDensity.nGridFreq
    error = onb.onbMatrix*(squeeze(specDensity.ONB(k,:,:)) - squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    m(k) = max(abs(error(:)));

    specTrue = onb.onbMatrix*(squeeze(specDensity_true.ONB(k,:,:)))*onb.onbMatrix';
    m2(k) = max(abs(specTrue(:)));
end
specDensity_error.mse_maxmax = max(m);
specDensity_error.rmse_maxmax = max(m)/max(m2);


end