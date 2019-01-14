function [ expWu] = expSE3(Wu )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find exponential map from se(3) to SE(3).
% dl= (w,u) in se(3) which us R^6.
% exp(dl) = exp([ W u; 0 0])
%         = [exp(W) Vu ;0 1]
% W is skew symmetric matrix for w.
% where V= I +((1 − cos(th))/th^2)*W +((th − sin(th))/th^3)*W^2
% th=w'*w;
% A =sin(th)/th;
% B =(1-cos(th))/th^2;
% C =(1 − A)/th^2;
% exp(W) = I + AW + BW^2
%       V= I + BW + CW^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: W =[W1;W2;W3;...], Wi is skew symmetric matrix and size (3,3).
%        u =[u1;u2;u3;...], ti is translation vectors and of size (3,1).
% Wui=[] 
% Output: exp(W,u)=[exp(W] | Vu ]
%                  [  0    |  1 ]
% expWu is (4*n,4).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[W,u]=SE3toRt(Wu);
nf=size(W,1)/3;

matW=vec2mat(W,9);
w=matW(:,[8,3,4]);
th=sum(w.^2,2).^(.5); % also norm_vecW
I=[ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1)];
WW=vec2mat(CompositionOfTransformation(W,W),9);

A=sin(th)./th;
B =(1-cos(th))./(th.^2);
C =(1-A)./(th.^2);

% expW=expS03(W);
expW= I + repmat(A,1,9).*matW +repmat(B,1,9).*WW;
V= I + repmat(B,1,9).*matW +repmat(C,1,9).*WW;

%%% handle numerical error
W_is_0=double(th<.0000001);
idx=find(W_is_0);
V(idx,:)=I(idx,:); 
V=vec2mat(V,3);
u=vec2mat(u,3);
Vu=ApplyTransformation(V,u);

%%
expW(idx,:)=I(idx,:); 
expW=vec2mat(expW,3);
expWu=RttoSE3(expW,Vu);
end
