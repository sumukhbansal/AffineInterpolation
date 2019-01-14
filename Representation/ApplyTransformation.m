function [ TX ] = ApplyTransformation(T,X)
% Apply Transformation Ti on Xi for i=1:N, where T: (3xN,3) and X: (N,3)
% TX=(T*X') can be written as follows for [3,3] matrices.
% TX=[T1;T2].*[X1;X1;X1;X2;X2;X2] , T1:(3,3) & X1 :(1,3)

[~,ns]=size(T);

switch ns
    %for transformation of size [3,3]
    case 3  
        TX=sum(T.*repelem(X,3,1),2);
        TX=vec2mat(TX,3);

    % for transformation of size [4,4]
    case 4
        X=[X, ones(size(X,1),1)];
        TX=sum(T.*repelem(X,4,1),2);
        TX=vec2mat(TX,4);
        TX=TX(:,1:3);
    otherwise
        disp('Wrong size of transformation');
     
end
end

