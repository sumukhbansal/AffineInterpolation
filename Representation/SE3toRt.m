function [ R,t] = SE3toRt( SE )
% convert SE3 to R,t;
% SE =[SE1;SE2;...], SEi is [4,4]
%  R =[R1;R2;R3,...], Ri is [3 3]
%  t =[t1;t2;t3,...], ti is [3,1]
SE=vec2mat(SE,16);
R=vec2mat(SE(:,[1 2 3 5 6 7 9 10 11]),3);
t=vec2mat(SE(:,[4 8 12]),1);

end

