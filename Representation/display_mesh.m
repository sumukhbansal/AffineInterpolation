function display_mesh( V,F,meshtype,cl)

lc=strcat('.',cl);
switch meshtype
    case 'mesh'
        p=patch('Vertices',V,'faces',F);
        set(p,'FaceColor',cl,'EdgeColor','black')
        axis image;
        %         axis square;
        %         camlight;
               axis off;
        %            grid on;
        cameratoolbar('setmode', 'orbit');
        cameratoolbar('setcoordsys', 'y');
        cameratoolbar('show');
        
    case 'points'
        plot3(V(:,1),V(:,2),V(:,3),lc);
        
end

