function [SE] =RttoSE3( R,t)
% convert  R,t to SE3;
% SE =[SE1;SE2;...], SEi is [4,4]
%  R =[R1;R2;R3,...], Ri is [3 3]
%  t =[t1;t2;t3,...], ti is [3,1]
n=size(R,1)/3;
if size(t,2)~=1
    t=vec2mat(t,1);
end
SE=vec2mat([R,t],12);
SE=vec2mat([SE zeros(n,3) ones(n,1)],4);

end