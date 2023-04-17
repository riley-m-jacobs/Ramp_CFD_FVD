function meshplot(u,E,V,colorlabel,heading,zoom,nodat)

% Creates a contour plot of data on the triangular control volumes. Need a
% mesh that has defined the edges (E) and V (vertices) of each triangle.
%
% --inputs--
% u          = data that is to be plotted... M or P
% E          = nodes of each of the triangles in the mesh
% V          = vertices of the mesh, V(:,1) = x and V(:,2) = y
% colorlabel = string that assigns the colorbar scale --> downloadthe file!
% heading    = basically the title of the plot
% zoom       = used to zoom into the 'origin', attempted for the wedge but not currently used 
% nodat      = 'Nodat' has been used to plot just the mesh, can be defined
%               to make it easier to visualize
%         ...


hold on

set(gca,'FontSize',12)

% Unstretch the axis
ax0 = [min(V(:,1)) max(V(:,1)) min(V(:,2)) max(V(:,2))]; 
axis(ax0);

% Set zoomed axis 
ax = axis./(zoom);
axis(ax);

cbar = colorbar;
cbar.Label.String = colorlabel;
xlabel('x Position'); 
ylabel('y Position'); 
title(heading,'Interpreter','none');

% Get approx plotting limits for zoom axis
limx_neg = ax(1)*1.75;
limx_pos = ax(2)*1.75;
limy_neg = ax(3)*1.75;
limy_pos = ax(4)*1.75;

% Finite volume of triangles:
for i = (1:length(E(:,1)))
    
    tri = zeros(3,2); 
    
    p1 = E(i,1); % Data of u to E
    p2 = E(i,2);
    p3 = E(i,3); 
     
    tri(1,1) = V(p1,1); % Assign x coordinates
    tri(2,1) = V(p2,1);
    tri(3,1) = V(p3,1);
 
    tri(1,2) = V(p1,2); % Assign y coordinates
    tri(2,2) = V(p2,2);
    tri(3,2) = V(p3,2);
    
    if min(tri(:,1)) >= limx_neg && max(tri(:,1)) <= limx_pos && min(tri(:,2)) >= limy_neg && max(tri(:,2)) <= limy_pos
        if nodat
            fill(tri(:,1),tri(:,2),[1 1 1]);
        else
            fill(tri(:,1),tri(:,2),u(i),'LineStyle','none');
        end
    end
end

hold off

end