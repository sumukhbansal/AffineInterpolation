function [iSE]=invSE3(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inverse of SE3.
%  Input :  SE=[SE1;SE2;SE3;...] SEj is in SE(3),
%  Output: iSE=[iSE1;iSE2;iSE3;...], iSEj is in SE(3), 
%  SE = [R t;0 1]
%  iSE= [R' -R'*t;0 1];
% R=[R1;R2;R3];
% t=[t1;t2;t3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin==1)
    SE=cell2mat(varargin(1));
    [R,t]=SE3toRt(SE);
else
    R=cell2mat(varargin(1));
    t=cell2mat(varargin(2));
end


if size(t,2)==1
    t=vec2mat(t,3);
end
% Transpose is inverse for rotation matrix.
iR=transposeMatrix(R);
it=-1*ApplyTransformation(iR,t);
iSE=RttoSE3(iR,it);
end