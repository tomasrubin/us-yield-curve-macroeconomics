function onb = fPrepareONB(nGridSpace, nGridSpace_spacing, nKnots)
% this function creates onb struct
% nGridSpace ... how many points are on spation grid
% nKnots ... how many knots should the B-spline basis have

% spacial grid for visualitation
onb.nGridSpace = nGridSpace;
onb.gridSpace = linspace(0,1,onb.nGridSpace);

% prepare the B-spline onb.basis
% addpath('fdaM')
onb.nKnots = nKnots;
onb.knots = linspace(0,1,onb.nKnots)';
onb.nOrder = 3;
onb.nBasis = length(onb.knots) + onb.nOrder - 2; % should be equal to onb.nGridSmoother 
onb.basis = create_bspline_basis([0,1],onb.nBasis,onb.nOrder,onb.knots);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grid for the smoothing procedure
onb.gs_indx = 1:nGridSpace_spacing:onb.nGridSpace; % XXX here it should match !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
onb.gridSmoother = onb.gridSpace(onb.gs_indx);
onb.nGridSmoother = length(onb.gridSmoother);  % should be equal to onb.nBasis

if onb.nGridSmoother ~= onb.nBasis
    fprintf('Error: onb.nGridSmoother ~= onb.nBasis\n')
    return;
end

% the "censor" matrix from space to smoother grid
onb.H_smoother2space = zeros(onb.nGridSmoother,onb.nGridSpace);
for ii = 1:onb.nGridSmoother
    onb.H_smoother2space(ii, onb.gs_indx(ii)) = 1;
end
onb.basisMatrix = full(eval_basis(onb.basis, onb.gridSpace, 0)); % in columns there are onb.basis functions, onb.basisMatrix(:,i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orthogonalize the B-spline onb.basis

[Q,R] = qr(onb.basisMatrix);

onb.transform_O2B = R(1:onb.nBasis, 1:onb.nBasis); % if you multiply from right by this matrix you go from the new onb.basis (orthonormal) to the old one
onb.transform_O2B = onb.transform_O2B ./ sqrt(length(onb.gridSpace)); % because the scalar product is the integral over [0,1]
onb.transform_B2O = inv(onb.transform_O2B); % multiplying by this goes from the B-spline onb.basis to the orthonormal
onb.onbMatrix = onb.basisMatrix * onb.transform_B2O; % in columns there are onb.basis functions, onb.basisMatrix(:,i)

% the matrix for interpolation of the onb.gridSmoother data but express it in the ONB
mmmmmm = onb.H_smoother2space*onb.onbMatrix;
onb.smoothingMatrix_gridSmoother2ONB = inv(mmmmmm'*mmmmmm)*mmmmmm';



end