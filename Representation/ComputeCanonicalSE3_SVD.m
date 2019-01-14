function SE=ComputeCanonicalSE3_SVD(PV,F)
% Compute the canonical SE3 

P1=vec2mat(PV(F(:,1),:),1);
P2=vec2mat(PV(F(:,2),:),1);
P3=vec2mat(PV(F(:,3),:),1);

PEdge1=PV(F(:,2),:)-PV(F(:,1),:); %V2-V1
PEdge2=PV(F(:,3),:)-PV(F(:,1),:); %V3-V1


% unit edge vectors
E1=NormalizeVector(PEdge1);
E2=NormalizeVector(PEdge2);


E1_x_E2=NormalizeVector(cross(E1,E2,2));

ABt=[vec2mat(E1,1),zeros(size(F,1)*3,1),vec2mat(E1_x_E2,1)];

Rcanonical=zeros(size(ABt));

for i=1:size(F,1)
    [u,~,v]=svd(ABt(3*(i-1)+1:3*(i-1)+3,:));
    d=[1 0 0;0 1 0;0 0 det(v*u')];
    Rcanonical(3*(i-1)+1:3*(i-1)+3,:)=v*d*u';
end

T=-ApplyTransformation(Rcanonical,vec2mat(P1,3));
SE=RttoSE3(Rcanonical,T);

end

