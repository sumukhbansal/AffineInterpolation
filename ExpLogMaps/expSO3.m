function [ expW] = expSO3( W )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find exponential map from so(3) to SO(3).
% exp : W to R, W=[w]+ in so(3) and R is in SO(3).
% W is skew symmetric matrix for w.
% th=w'*w;
% A =sin(th)/th;
% B =(1-cos(th))/th^2;
% C =(1 − A)/th^2;
% expSO3(W) = I + AW + BW^2
% References: 
% https://pixhawk.org/_media/dev/know-how/jlblanco2010geometry3d_techrep.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matW=vec2mat(W,9);
w=matW(:,[8,3,4]);
nf=size(matW,1);
theta=sum(w.^2,2).^(.5);
% theta=sum(w.^2,2); % theta
I=[ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1)];
WW=vec2mat(CompositionOfTransformation(W,W),9);

A=sin(theta)./theta;
B =(1-cos(theta))./(theta.^2);
expW= I + repmat(A,1,9).*matW +repmat(B,1,9).*WW;

%%% handle numerical error
W_is_0=double(theta<.0000001);
idx=find(W_is_0);
expW(idx,:)=I(idx,:); 
expW=vec2mat(expW,3);



end
