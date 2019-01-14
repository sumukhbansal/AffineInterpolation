function [TSEp ] = logSE3(expSEp )
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Log map from SE(3) to se(3).
% dl= (w,u) in se(3) which us R^6.
% exp(dl) = exp([ W u; 0 0])
%         = [exp(W) Vu ;0 1]
% W is skew symmetric matrix for w.
% where V= I +((1 − cos(th))/th^2)*W +((th − sin(th))/th^3)*W^2
% th=w'*w;
% A =sin(th)/th;
% B =(1-cos(th))/th^2;
% C =(1 − A)/th^2;
% W = logS03(expW)
% invV= I - W/2 +(1/th^2)*(1 −A/2B)W^2
% u= invV* Vu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: W =[expW1;expW2;expW3;...], expWi is skew symmetric matrix and size (3,3).
%        u =[Vu1;Vu2;Vu3;...], Vui is translation vectors and of size (3,1).
% Output: log(expW,Vu)=[logSO3(expW) | V^-1 *Vu ]
%                      [  0          |  0       ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[expW,Vu]=SE3toRt(expSEp);
W=logSO3( expW );

matW=vec2mat(W,9);
w=matW(:,[8,3,4]);
nf=size(w,1);
th=sum(w.^2,2).^(.5); % also norm_vecW
I=[ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1) zeros(nf,1) zeros(nf,1) zeros(nf,1) ones(nf,1)];
WW=vec2mat(CompositionOfTransformation(W,W),9);

A=sin(th)./th;
B =(1-cos(th))./(th.^2);

invV= I - .5*matW + repmat((1./th.^2),1,9).*(1-repmat(A./(2*B),1,9)).*WW;

% check numerical error
W_is_0=double(th<.0000001);
idx=find(W_is_0);
invV(idx,:)=I(idx,:); 
invV=vec2mat(invV,3);
if size(Vu,2)==1
    Vu=vec2mat(Vu,3);
end
u=vec2mat(ApplyTransformation(invV,Vu),1);
TSEp=vec2mat([vec2mat([W,u],12),zeros(nf,4)],4);

end

