function expsym=expSYM3(Y)

[~,ev]=eig(Y);
l=sort(diag(ev));
l3=l(1);l2=l(2);l1=l(3);
Z=Y-l2*eye(3);
l11=l1-l2;l22=0;l33=l3-l2;
b=1-( (l11*l22*(e2(l11)-e2(l33)))/(l11-l33));
c=.5+( (l11*(2*e2(l11)-1) - l33*(2*e2(l33)-1))/(2*(l11-l33)));
expsym=exp(l2)*(eye(3)+b*Z+c*Z*Z);
end

function e2x= e2(x)
e2x=(exp(x)-1-x)/(x*x);
end