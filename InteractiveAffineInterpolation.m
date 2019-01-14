% Test Affine Interpolation
clear all;
close all;
clc;

%% Paths for dependencies
addpath(genpath('ExpLogMaps'));
addpath(genpath('Representation'));


%% Take input Affine transformation from the user
d = input('Enter 3x1 translation part of the Affine transformation.\n');
d=vec2mat(d,1);
L = input('Enter 3x3 Linear part of the Affine transformation.\n');
while det(L)<.01
    disp("Invalid orientation preserving transformation or Determinant is very small.\n");
    L = input('Enter 3x3 Linear part of an orientation preserving affine transformation.\n');
end

clc;
disp('Enter the object on which the interpolation is to be shown.')
type = input('Enter 1 for tetrahedron and 2 for cube.\n');
clc;
timestep = input('Enter the time step for the interpolation like [.25 .3 .8] or 0:.1:1 \n');
clc;
disp('Affine transformation T =');disp([L,d; 0 0 0 1])
if type==1
    disp('Interpolation is to be shown on tetrahedron.');
else
    disp('Interpolation is to be shown on cube.');
end
disp('timestep =');
disp(timestep);
disp('Press a key for the interpolation');
pause();
%% source and target tetrahedron
I=[0 0 0; 1 0 0; 1 1 0;0 0 1];
P1=[0 0 0;1 0 0;1 1 0;0 0 1];
P2=(L*P1')'+repmat(d',4,1);

%% compute the Lie group representation
F=[1 2 3 4];
[Rp, Ap, Sp]=Compute_EAS(P1,I,F);
[Rq, Aq, Sq]=Compute_EAS(P2,I,F);
[TR, TA, TS]=logEAS(Rp,Ap,Sp,Rq,Aq,Sq);

Yc=zeros(length(timestep),12);
k=1;
for i=timestep
    TRn=i*TR; TAn=i*TA; TSn=i*TS;
    [Rc, Ac, Sc]=expEAS(Rp,Ap,Sp,TRn,TAn,TSn);
    Tc=ComputeT_EAS(Rc,Ac,Sc,I);
    Yc(k,:)= vec2mat(ShapeReconstruction_EAS(Tc, I, F),12);
    k=k+1;
end

switch type
    case 1 %Tetrahedron
        F=[1 2 3; 2 3 4; 3 1 4; 1 2 4];
        for i=1:length(timestep)
            display_mesh(vec2mat(Yc(i,:),3),F,'mesh','y'); hold on;
        end
        display_mesh(P1,F,'mesh','g'); hold on;
        display_mesh(P2,F,'mesh','g'); hold on;
    case 2 %Cube
        
        Vi=[-.5,-.5,.5;-.5,.5,.5;.5,.5,.5;.5,-.5,.5;-.5,-.5,-.5;-.5,.5,-.5;.5,.5,-.5;.5,-.5,-.5];
        F=[2 1 4;4 3 2;3 4 8;8 7 3;4 1 5;5 8 4;7 6 2;2 3 7;5 6 7;7 8 5;6 5 1;1 2 6];
        
        for i=1:length(timestep)
            [Af]=ComputeAffinity(P1,vec2mat(Yc(i,:),3));Vc= (Af.L*Vi')'+repmat(Af.O',8,1);
            display_mesh(Vc,F,'mesh','y'); hold on;
        end
        [Af]=ComputeAffinity(P1,P2);Vf= (Af.L*Vi')'+repmat(Af.O',8,1);
        display_mesh(Vf,F,'mesh','g'); hold on;
        display_mesh(Vi,F,'mesh','g'); hold on;
    otherwise
        disp('Case not valid');
end

