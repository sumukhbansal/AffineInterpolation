function T = ComputeT_EAS( SEq,Aq,Sq,IV )
%% compute matrix T from (SEq Aq Sq) 
% IV is the template simplex.
n=size(SEq,2);
nf=size(Sq,1);
Aq=vec2mat([vec2mat([Aq,zeros(3*nf,1)],12), zeros(nf,3) ones(nf,1)],4);
SEx=repmat(ComputeCanonicalSE3_SVD(IV,[1 2 3]),nf,1);
Sq=vec2mat(kron(Sq,eye(4)),16);
Sq(:,16)=1;
Sq=vec2mat(Sq,4);
SSEx=CompositionOfTransformation(Sq,SEx);
ASSEx=CompositionOfTransformation(Aq,SSEx);
iSEx=invSE3(SEx);
iSExASSEx=CompositionOfTransformation(iSEx,ASSEx);
T=CompositionOfTransformation(SEq,iSExASSEx);
end
