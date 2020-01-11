function y_out = matlab_viz(A, cmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2019, Orit Peleg, orit.peleg@colorado.edu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dbstop if error
set(0,'Defaultlinelinewidth',3.5, 'DefaultlineMarkerSize',12,...
    'DefaultTextFontSize',5, 'DefaultAxesFontSize',12);
global scrsz
scrsz = get(0,'ScreenSize');

set(0, 'DefaultFigureRenderer', 'OpenGL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
print_3D_fancy(A, cmax); 
y_out = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function y=print_3D_fancy(X, cmax)
 
r = 0.5;
x = X;
[X,Y,Z]=cylinder(1,4);

cmap = colormap(parula(cmax));
for i=1:length(x(:,1))
    h = surf(X*r+x(i,1),Y*r+x(i,2),Z*r*2+x(i,3));
    z_i=int8(x(i,3))+1;
    if z_i>cmax
        z_i = cmax;
    end
    set(h, 'FaceColor',cmap(z_i,1:3));
    hold on;
end

xlabel('x'); ylabel('y'),zlabel('z');
 
material shiny
lighting gouraud
alpha(0.99);
axis equal;
xlim([-10 10]); ylim([-10 10]); zlim([0,cmax])

 
view(10,30);
camlight(10,5);
 
y = 0;
end
 
