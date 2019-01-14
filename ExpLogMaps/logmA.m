function lA= logmA(A)
a1=A(1,2);
a2=A(2,2);
a3=A(1,3);
a4=A(2,3);
a5=A(3,3);
type=0;
if (abs(a2-1)<=.00001 && abs(a5-1)>0.00001)
    type=1;
end
if (abs(a5-1)<=.00001 && abs(a2-1)>0.00001)
    type=2;
end
if (abs(a2-a5)<.00001 && abs(a5-1)<=0.00001)
    type=3;
end
if (abs(a2-a5)<.00001 && abs(a5-1)>.00001)
    type=4;
end
switch type
    case 1
        lA= [ 0,        a1,                   (log(a5)*(a1*a4 - a3 + a3*a5))/(a5^2 - 2*a5 + 1) - (a1*a4)/(a5 - 1);
             0,         0,                    (a4*log(a5))/(a5 - 1);
             0,         0,                     log(a5)                                                           ];

    case 2
        lA= [ 0,        (a1*log(a2))/(a2-1),        (a1*a4*log(a2))/(a2^2 - 2*a2 + 1) - (a3 + a1*a4 - a2*a3)/(a2 - 1);
             0,          log(a2),                   (a4*log(a2))/(a2 - 1);
             0,          0,                         0                                                                ];

    case 3
        lA= [ 0,              a1,                       a3- (a1*a4)/2 ;
             0,               0,                        a4;
             0,               0,                        0                       ];

    case 4
        lA= [ 0, (a1*log(a2))/(a2-1),  (a1*a4*(1/a2 + (a3*log(a2))/(a1*a4)))/(a2 - 1) - (a1*a4*log(a2))/(a2^2 - 2*a2 + 1);
             0,               log(a2),                                  (a4/a2);
             0,                     0,                                                                        log(a2)];

    otherwise
       lA= [ 0, (a1*log(a2))/(a2-1), (a1*a4*log(a2))/((a2-a5)*(a2-1)) - (log(a5)*(a1*a4-a2*a3+a3*a5))/((a2-a5)*(a5-1));
             0,               log(a2),                                  (a4*log(a5/a2))/(a5-a2);
             0,                     0,                                                                        log(a5)];
 
end


end