function rmse = fKrigingZ_error( response_est, response, var_z )

rmse = mean( (response-response_est).^2 )/var_z;

end