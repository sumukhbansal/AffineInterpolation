function [ TSEpq,TApq,TSpq ] = logEAS( SEp,Ap,Sp,SEq,Aq,Sq )
TSEpq=logSE3(CompositionOfTransformation(SEq,invSE3(SEp)));
TSpq=logS(Sq./Sp); 
TApq=logmA(CompositionOfTransformation(inv(Ap),Aq));  
end
