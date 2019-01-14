function [ SEq,Aq,Sq ] = expEAS( SEp,Ap,Sp,TSE,TA,TS )
SEq=CompositionOfTransformation(expSE3(TSE),SEp);
Aq=CompositionOfTransformation(Ap,expm(TA));
Sq=Sp.*expS(TS);
end
