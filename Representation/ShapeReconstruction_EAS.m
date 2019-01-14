function [ Y ] = ShapeReconstruction_EAS( T, IV,F )

[n]=size(T,2);
nf=size(F,1);
nv=max(max(F));

PEdge0=[repmat(IV(1,:),nf,1),ones(nf,1)];
PEdge1=[repmat(IV(2,:),nf,1),ones(nf,1)];
PEdge2=[repmat(IV(3,:),nf,1),ones(nf,1)];
PEdge3=[repmat(IV(4,:),nf,1),ones(nf,1)];

E1=F(:,1);
E2=F(:,2);
E3=F(:,3);
E4=F(:,4);
I=[1:nf]';
K=double(ones(nf,1));
A1 = sparse(I,E1,K,nf,nv);
A2 = sparse(I,E2,K,nf,nv);
A3 = sparse(I,E3,K,nf,nv);
A4 = sparse(I,E4,K,nf,nv);
AE1=kron(A1,eye(n));
AE2=kron(A2,eye(n));
AE3=kron(A3,eye(n));
AE4=kron(A4,eye(n));

% break down for low memory
nf1=4000;
Diag_T=[];
for st=1:nf1:size(T,1)
    nf1=min(size(T,1)-st+1,nf1);
    T1=T(st:st+nf1-1,:);
    temp1=sparse(repmat(T1,1,nf1/4));
    temp2=kron(sparse(eye(nf1/4)),ones(n));
    temp_Diag_T=temp1.*temp2;
    Diag_T=sparse(blkdiag(Diag_T,temp_Diag_T));
end
AE=sparse([AE1;AE2;AE3;AE4]);
b=sparse([Diag_T*vec2mat(PEdge0,1);Diag_T*vec2mat(PEdge1,1);Diag_T*vec2mat(PEdge2,1);Diag_T*vec2mat(PEdge3,1)]);
% tic;
Y=full(AE\b);
% tm=toc;
% fprintf('Computing Least square. Time taken: %d \n',tm);
Y=vec2mat(Y,n);
Y=Y(:,1:3);

end

