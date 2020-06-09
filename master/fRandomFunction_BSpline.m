function f_inSpace = fRandomFunction_BSpline(onb,nGridSpace,smoothness,norm_target)

% This function generates a random smooth function in the spatial domain
% nGridSpace ... the number of uniform grid points
% smoothness ... integer - a parameter how much smooth, more = smoother
% norm_target .. positive real - what the L^2([0,1]) norm of the random operator should be
%   norm(m_inSpace) / sqrt(nGridSpace) ... is set to be 'norm_target'

% example call: fRandomFunction_Fourier(241,1,1)


% create Fourier basis
addpath('fdaM')
nBasis = 101;
basis = create_fourier_basis([0,1],nBasis);
basisMatrix = full(eval_basis(basis, linspace(0,1,nGridSpace), 0));


f_Fourier = nan(nBasis,1);
for ii = 1:nBasis    
    f_Fourier(ii) = randn(1) * exp( -0.5*smoothness*(ii) );
end

% plot(basisMatrix*f_Fourier)


f_inSpace = basisMatrix*f_Fourier;

% make sure that the norm of the operator is as wanted
f_inSpace = f_inSpace / norm(f_inSpace) * norm_target * sqrt(nGridSpace);

end