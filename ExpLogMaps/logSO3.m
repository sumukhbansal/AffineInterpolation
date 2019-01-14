function [ logR ] = logSO3( R )
% find exponential map of an element in SO3.
% ln(R) = (θ/(2 sin θ)) (R − R')
% cos θ =(tr(R) − 1)/2
% ω = [ln(R)]▽
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matR=vec2mat(R,9);
traceR=sum(matR(:,[1 5 9]),2);
nf=size(matR,1);
theta=acos(((traceR-1)/2));

%check numerical error
theta_is_0=double(abs(theta)<.0000001);
idx=find(theta_is_0);
A=sin(theta)./theta;
A(idx,:)=1; 


logR=repmat(1./(2*A),1,9).*vec2mat(R-transposeMatrix(R),9);

R_is_I=sum(abs(logR),2)<.0000001;
idx1=find(R_is_I);
W=.0001*[zeros(nf,1) -ones(nf,1) ones(nf,1) ones(nf,1) zeros(nf,1) -ones(nf,1) -ones(nf,1) ones(nf,1) zeros(nf,1)];
logR(idx1,:)=W(idx1,:); 
logR=vec2mat(logR,3);

end
