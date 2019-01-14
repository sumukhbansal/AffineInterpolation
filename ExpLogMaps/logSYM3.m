function lgsym= logSYM3(S)
tAA=S*S;
[~,ev]=eig(tAA);
l=sort(diag(ev));
l3=l(1);l2=l(2);l1=l(3);
Z=tAA/l2;
l11=l1/l2; l22=1; l33=l3/l2;
a=-1+((l33*L2(l11)-l11*L2(l33))/(l11-l33));
c=(L2(l11)-L2(l33))/(l11-l33);
lgsym=((a+log(l2))*eye(3)-(a+c)*Z+c*Z*Z)/2;
end

function l2x=L2(x)
l2x=(log(x)-(x-1))/(x-1);
end