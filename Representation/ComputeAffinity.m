%%% Compute affinity from tetrahedron pairs or triangle pairs
function A=ComputeAffinity(P,Q)
n=size(P,1)-1;
ok=0;
aux=[P,ones(n+1,1)]';
if (abs(det(aux))<1.0e-10)
    L = eye(n,n); O=zeros(n,1); %iniId();
    ok=1;
end
r=Q';
res=r*(inv(aux));
L=res(1:n,1:n);
O=res(1:n,n+1);
if (abs(det(L))<1.0e-10)
    ok=2;
end
A.L=L;
A.ok=ok;
A.O=O;
end
