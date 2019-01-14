function RT=transposeMatrix(R)

n=size(R,2);
switch n
    case 3
        VR=vec2mat(R,9);
        VR=VR(:,[1 4 7 2 5 8 3 6 9]);
        RT=vec2mat(VR,3);
    case 4
        VR=vec2mat(R,16);
        VR=VR(:,[1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16]);
        RT=vec2mat(VR,4);
    case 1
        RT=R;
end
end