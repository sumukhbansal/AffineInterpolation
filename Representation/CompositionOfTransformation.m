function [ T ] = CompositionOfTransformation(T2,T1)
% Combine two arrays of transformations.
% Ti=[Ri1;...;Rin]
% T=[Rj1*Ri1 ;...; Rjn*Rin], j=2 i=1;

if (size(T1,2)~=size(T2,2))
    disp('mismatch in Dimentions');
end
n=size(T1,2);
switch n
    case 3
        % Mat version of each column of T1;
        T11=vec2mat(T1(:,1),3);
        T12=vec2mat(T1(:,2),3);
        T13=vec2mat(T1(:,3),3);
        
        T21=sum(T2.*repelem(T11,3,1),2);
        T22=sum(T2.*repelem(T12,3,1),2);
        T23=sum(T2.*repelem(T13,3,1),2);
        
        T=[T21 T22 T23];
    case 4
        % Mat version of each column of T1;
        T11=vec2mat(T1(:,1),4);
        T12=vec2mat(T1(:,2),4);
        T13=vec2mat(T1(:,3),4);
        T14=vec2mat(T1(:,4),4);
        
        T21=sum(T2.*repelem(T11,4,1),2);
        T22=sum(T2.*repelem(T12,4,1),2);
        T23=sum(T2.*repelem(T13,4,1),2);
        T24=sum(T2.*repelem(T14,4,1),2);
        
        T=[T21 T22 T23 T24];
        
    case 1
        T=T2.*T1;
    otherwise
        disp('Composition not defined');
end





end