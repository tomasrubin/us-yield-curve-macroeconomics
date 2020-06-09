function m_inSpace = fRandomOperator(nGridSpace,symetry,smoothness,norm_target)

% This function generates a random smooth operator in the spatial domain
% nGridSpace ... the number of uniform grid points
% symetry ...... 1 if the operator should be symetric & positve semi-definite, 0 if not
% smoothness ... integer - a parameter how much smooth, more = smoother
% norm_target .. positive real - what the Hilbert-Schmidt (aka Frobenius) norm of the random operator should be
%   norm(m_inSpace) / nGridSpace ... is set to be 'norm_target'

if symetry % symetric operator
    m_inSpace = randn(nGridSpace);
    m_inSpace = m_inSpace*m_inSpace';
    smoothMat = fSmoothMat( nGridSpace, smoothness);
    m_inSpace = smoothMat * m_inSpace * smoothMat';
    m_inSpace = (m_inSpace+m_inSpace')/2;    
else % nonsymetric
    m_inSpace = randn(nGridSpace);
    smoothMat = fSmoothMat( nGridSpace, smoothness);
    m_inSpace = smoothMat * m_inSpace * smoothMat';
end

% make sure that the norm of the operator is as wanted
m_inSpace = m_inSpace / norm(m_inSpace) * norm_target * nGridSpace;

end