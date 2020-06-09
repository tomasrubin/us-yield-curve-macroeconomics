function m_inSpace = fRandomOperator_Fourier(nGridSpace,symetry,smoothness,exponent,norm_target)

% This function generates a random smooth operator in the spatial domain
% nGridSpace ... the number of uniform grid points
% symetry ...... 1 if the operator should be symetric & positve semi-definite, 0 if not
% smoothness ... integer - a parameter how much smooth, more = smoother
% exponent ... positive real - what is the exponent of (ii) or (ii+jj)
% norm_target .. positive real - what the Hilbert-Schmidt (aka Frobenius) norm of the random operator should be
%   norm(m_inSpace) / nGridSpace ... is set to be 'norm_target'

% example call: fRandomOperator_Fourier(241,1,1,1)

% create Fourier basis
addpath('fdaM')
nBasis = 101;
basis = create_fourier_basis([0,1],nBasis);
basisMatrix = full(eval_basis(basis, linspace(0,1,nGridSpace), 0));

m_Fourier = nan(nBasis);
if symetry % symetric
    
    % random operator in Fourier basis
    for ii = 1:nBasis
        for jj = 1:nBasis
            m_Fourier(ii,jj) = randn(1) * exp( -0.5*smoothness*(ii)^exponent );
        end
    end
    
    % convert into spatial domain
    m_inSpace = basisMatrix*m_Fourier*basisMatrix';    
    m_inSpace = m_inSpace*m_inSpace';    
    m_inSpace = (m_inSpace+m_inSpace')/2;    
    
else % asymetric
    
    % random operator in Fourier basis
    for ii = 1:nBasis
        for jj = 1:nBasis
            m_Fourier(ii,jj) = randn(1) * exp( -0.5*smoothness*(ii+jj)^exponent );
        end
    end    
    % convert into spatial domain
    m_inSpace = basisMatrix*m_Fourier*basisMatrix';
end

% make sure that the norm of the operator is as wanted
m_inSpace = m_inSpace / norm(m_inSpace) * norm_target * nGridSpace;
% surf(m_inSpace)
% eig(m_inSpace)

end