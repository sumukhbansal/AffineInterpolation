function [Y] = NormalizeVector(X)
% normalize array of vectors row wise. To normalize edges and normals
% X : nx3

normX=sqrt(sum(X.*X,2));

% solve the devide by zero error
normX=normX+double(normX<eps);

Y=X./repmat(normX,1,size(X,2));

end

