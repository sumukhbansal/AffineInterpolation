function [SE,A,S]=Compute_EAS(QV,IV,F1)
% QV is the tetrahedron to be represnted wrt template tetrahedra IV.
% QV and IV both contains 4 vertices of the tetrahedron.
% F represent the seqence of the vertices which decide the first face. 

F=F1(:,1:3);
nf=size(F,1);

P1=repmat(IV(1,:),nf,1);
P2=repmat(IV(2,:),nf,1);
P3=repmat(IV(3,:),nf,1);
P4=repmat(IV(4,:),nf,1);


%% Compute canonical SE matrix for both poses
SEp=repmat(ComputeCanonicalSE3_SVD(IV,[1 2 3]),nf,1);
SEq=ComputeCanonicalSE3_SVD(QV,F);

P11=ApplyTransformation(SEp,P1);
P22=ApplyTransformation(SEp,P2);
P33=ApplyTransformation(SEp,P3);
P44=ApplyTransformation(SEp,P4);

Q1=QV(F1(:,1),:);
Q2=QV(F1(:,2),:);
Q3=QV(F1(:,3),:);
Q4=QV(F1(:,4),:);


Q11=ApplyTransformation(SEq,Q1);
Q22=ApplyTransformation(SEq,Q2);
Q33=ApplyTransformation(SEq,Q3);
Q44=ApplyTransformation(SEq,Q4);




S=(Q22(:,1)./P22(:,1));
 
%% Afine Transformation
P11=P11.*repmat(S,1,3);
P22=P22.*repmat(S,1,3);
P33=P33.*repmat(S,1,3);
P44=P44.*repmat(S,1,3);

% A =  [1 alpha1 alpha3 0;
%       0  alpha2 alpha4 0;
%       0   0     alpha5 0;
%       0   0       0    1;]
% such that A*[w1 w2']=[1 1 0 1;0 0 1 1]'

alpha5=Q44(:,3)./P44(:,3);
alpha2=Q33(:,2)./P33(:,2);
alpha1=(Q33(:,1)-P33(:,1))./(P33(:,2));
alpha4=(P33(:,2).*Q44(:,2)-P44(:,2).*Q33(:,2))./(P33(:,2).*P44(:,3));
alpha3=(P33(:,2).*Q44(:,1)-P44(:,2).*Q33(:,1) + P33(:,1).*P44(:,2)- P44(:,1).*P33(:,2))./(P33(:,2).*P44(:,3));


nf=size(F,1);
A=[ones(nf,1) alpha1 alpha3 zeros(nf,1) alpha2  alpha4 zeros(nf,1) zeros(nf,1) alpha5];
A=vec2mat(A,3);
 

iSEq=invSE3(SEq);
SE=CompositionOfTransformation(iSEq,SEp);

end
