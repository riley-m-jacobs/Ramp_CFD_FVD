function [extract_line] = bottom_plot(V, E, prop)

%
% This is the function for plotting along the bottom boundary
% --inputs--
% V          = vertices of the mesh, V(:,1) = x and V(:,2) = y
% E          = nodes of each of the triangles in the mesh
% prop       = property you want to investigate
%         ...
%
% --outputs--
% extract_line = extracts the values of the property you want along bottom
% boundary
%

ax0 = [min(V(:,1)) max(V(:,1)) min(V(:,2)) max(V(:,2))]; 

limx_neg = ax0(1)*1.75;
limx_pos = ax0(2)*1.75;
limy_neg = ax0(3)*1.75;
limy_pos = ax0(4)*1.75;

datamatrix = [];

for i = (1:length(E(:,1)))
    
    tri = zeros(3,2); 
    
    p1 = E(i,1); 
    p2 = E(i,2);
    p3 = E(i,3); 
     
    tri(1,1) = V(p1,1); 
    tri(2,1) = V(p2,1);
    tri(3,1) = V(p3,1);
 
    tri(1,2) = V(p1,2); 
    tri(2,2) = V(p2,2);
    tri(3,2) = V(p3,2);
   
    if min(tri(:,1)) >= limx_neg && max(tri(:,1)) <= limx_pos && min(tri(:,2)) >= limy_neg && max(tri(:,2)) <= limy_pos
        datamatrix = [datamatrix; tri(:,1)', tri(:,2)', prop(i)]; 
    end
end

extract_line = [datamatrix(1,1) datamatrix(1,4) datamatrix(1,7)];
j = 2;
for i = 1:(n_y-1)
    x_values = datamatrix(j,1);
    y_values= datamatrix(j,4);
    values= datamatrix(j,7);
    extract_line = [extract_line; x_values, y_values, values];
    j = j + 2;
end

end
