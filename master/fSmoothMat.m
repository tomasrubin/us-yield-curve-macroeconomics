function m = fSmoothMat(n,lambda)

m = zeros(n);
for i = 1:n
    for j = 1:n
        m(i,j) = exp(-(-i+j)^2/n^2/lambda) /20;
    end
end

% normalize to infinity norm 1
topE = max(eig(m));
m = m / topE;

end