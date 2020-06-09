function specCrossDensity_error = fSpecCrossDensityError( onb, specCrossDensity, specCrossDensity_true )

specCrossDensity_error = [];

%% L2 norm
specCrossDensity_error.mse_qv = 0;
divisor = 0;
for k = 1:specCrossDensity.nGridFreq
    m1 = onb.onbMatrix * ( specCrossDensity.ONB(k,:) - specCrossDensity_true.ONB(k,:) )';
    m2 = onb.onbMatrix * ( specCrossDensity_true.ONB(k,:) )';
    specCrossDensity_error.mse_qv = specCrossDensity_error.mse_qv + 2*pi/specCrossDensity.nGridFreq * sum( abs(m1).^2 );
    divisor = divisor + 2*pi/specCrossDensity.nGridFreq * sum( abs(m2).^2 );
end
specCrossDensity_error.rmse_qv = specCrossDensity_error.mse_qv / divisor;


%% intmax norm
specCrossDensity_error.mse_intmax = 0;
divisor = 0;
for k = 1:specCrossDensity.nGridFreq
    m1 = onb.onbMatrix * ( specCrossDensity.ONB(k,:) - specCrossDensity_true.ONB(k,:) )';
    m2 = onb.onbMatrix * ( specCrossDensity_true.ONB(k,:) )';
    specCrossDensity_error.mse_intmax = specCrossDensity_error.mse_qv + 2*pi/specCrossDensity.nGridFreq * max( abs(m1) );
    divisor = divisor + 2*pi/specCrossDensity.nGridFreq * max(abs(m2));
end
specCrossDensity_error.rmse_intmax = specCrossDensity_error.mse_intmax / divisor;




end