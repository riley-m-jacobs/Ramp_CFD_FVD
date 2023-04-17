function [E,V,B,n_y,bedges,inedges,x_ramp_start,x_length] = wedge_mesh(n_x,discretize,thetad)

%
% Generates a mesh for a wind tunnel with varying degrees of theta and length 
% depending on users choice. Can change the wind tunnel dimensions, 
% the discertization in the x-direction and along the ramp. N_y will match
% the number of grid points in the x so it will be an equal number of grid
% points in x and y. IN THE CODE: can change - x_min, x_max, y_min, y_max,
% x_ramp_start and x_length (of the ramp)
%
% --inputs--
% n_x        = number of grid points on the pre-shock and post-shock regions
% discretize = number of points + 2 that are up the ramp
% thetad     = theta of ramp in degrees
%         ...
%
% --outputs--
% E          = nodes of each of the triangles in the mesh
% V          = vertices of the mesh, V(:,1) = x and V(:,2) = y
% B          = boundary data with node correlations 
% n_y        = output a mesh that is n_y x n_y grid points
% bedges     = boundary edge data for the mesh as a matrix of the form [nA nB nx ny dl Eb] 
% inedges    = interior edge data for the mesh as a matrix of the form [nA nB nx ny dl Eb]
%
 
% Define the rectangular domain
x_min = 0;
x_max = 3;
y_min = 0;
y_max = 1;

% Define the number of nodes in the y directions to make equal numbered
% grid
n_y = (2*n_x) + (discretize-2);

% Create a grid of nodes
x = linspace(x_min, x_max, n_x);
y = linspace(y_min, y_max, n_y);
[X, Y] = meshgrid(x, y);

% Reshape the grid of nodes into a list of points
P = [X(:), Y(:)];

% Define the ramp geometry
x_length = 0.6;
%x_length = 0.31;
x_ramp_start = 0.5;
y_length = x_length*tand(thetad);

x_ramp = linspace(x_ramp_start, x_ramp_start + x_length, discretize);
y_ramp = linspace(0, y_length, discretize);

doublecheck = atand(y_length/x_length);
fprintf('Theta = %0.4f', doublecheck)

% Define the bounding box for the ramp
x_min_ramp = min(x_ramp);
x_max_ramp = max(x_ramp);
y_min_ramp = min(y_ramp);
y_max_ramp = max(y_ramp);
dx = x_max_ramp-x_min_ramp;

% Ramp quadrant
for i = 1:length(x_ramp)
    x_pt = x_ramp(i);
    y_slope(:,i) = linspace(y_ramp(i),y_max,n_y)';
end

for i = 1:size(y_slope,2)
    x_pt = x_ramp(i);
    x_slope(:,i) = x_pt*ones(1,size(y_slope,1))';
end

P_ramp = [x_slope(:) y_slope(:)];
P_ramp = P_ramp((n_y+1):end-n_y,:);

% Lower quadrant grid
%x_low = linspace(x_min,x_min_ramp,n_x*0.5);
x_low = linspace(x_min,x_min_ramp,n_x*0.25);

y_low = linspace(y_min,y_max,n_y);
[X_low, Y_low] = meshgrid(x_low, y_low);
P_low = [X_low(:),Y_low(:)];

% Upper quadrant grid
x_upp = linspace(x_max_ramp, x_max, n_x*1.75);
y_upp = linspace(y_max_ramp, y_max, n_y);
[X_upp, Y_upp] = meshgrid(x_upp, y_upp);
P_upp = [X_upp(:),Y_upp(:)];


P = [P_low; P_ramp; P_upp];

%scatter(P(:,1),P(:,2),'.','b')

x = reshape(P(:,1),[n_y,length(P(:,1))/n_y]);
y = reshape(P(:,2),[n_y,length(P(:,2))/n_y]);

% Streaching nodes towards the bottom profile
sf = 1.5; % set stretching factor
for i = 1:length(y(1,:))
    ymin = y(1,i);
    yn(:,i) = (y(:,i)-ymin)/(y_max-ymin);
end
%yn = (y-y_min)/(y_max-y_min);
ys = (1-exp(sf*yn))/(1-exp(sf));
for i = 1:length(y(1,:))
    ymin = y(1,i);
    y(:,i) = ys(:,i)*(y_max-ymin) + ymin;
end

mesh(x,y,ones(size(x))) % verify!

X = x;
Y = y;

x = X(1,:);
y = Y(:,1);

% Determine the total number of nodes in the mesh
numNodes = n_y^2;

% Initialize the Nodes array with zeros
Nodes = zeros(numNodes, 2);

% Assign the x and y coordinates to the Nodes array
for j = 1:n_y
    for i = 1:n_y
        c = n_y*(j-1) + i;
        Nodes(c, 1) = x(i);
        Nodes(c, 2) = y(j);
    end
end

% Determine the total number of triangles in the mesh
numTriangles = 2*(n_y-1)*(n_y-1);

% Initialize the triangles array with zeros
Triangles = zeros(numTriangles, 3);

% Assign the node indices to the triangles array
for j = 1:n_y-1
    for i = 1:n_y-1
        % Determine the indices of the four nodes that make up the rectangle
        node1 = n_y*(j-1) + i;
        node2 = n_y*(j-1) + i+1;
        node3 = n_y*j + i;
        node4 = n_y*j + i+1;
        
        % Assign the indices of the two triangles that make up the rectangle
        triangle1 = 2*(n_y-1)*(j-1) + 2*(i-1) + 1;
        triangle2 = triangle1 + 1;
        
        % Assign the node indices to the triangles array
        Triangles(triangle1, :) = [node1, node2, node3];
        Triangles(triangle2, :) = [node2, node4, node3];
    end
end

% Define the vertices of the triangles
E = Triangles;
V = zeros(length(P),2);

for i = 1:n_y
    for j = 1:n_y
        new = [X(i,j) Y(i,j)];
        V = [V; new];
    end
end

V = V(length(P)+1:end,:);

% Now define the different boundaries based on the nodes
B = cell(4,3);
lower = linspace(1,n_y,n_y);
top = linspace((n_y^2) - (n_y-1), (n_y^2), n_y);
right = linspace(n_y,n_y^2,n_y);
left = linspace(1,(n_y^2) - (n_y-1),n_y);

Bottom = zeros(n_y-1,2);
Top = zeros(n_y-1,2);

for i = 1:n_y-1
    node1 = lower(i);
    node2 = lower(i) + 1;
    Bottom(i,:) = [node1 node2];

    node3 = top(i);
    node4 = top(i) + 1;
    Top(i,:) = [node4 node3];

    node5 = left(i+1);
    node6 = left(i);
    Left(i,:) = [node5 node6];

    node7 = right(i);
    node8 = right(i+1);
    Right(i,:) = [node7 node8];

end

B{1,1} = 'Bottom';  B{1,2} = length(Bottom(:,1)); B{1,3} = Bottom;
B{2,1} = 'Right';   B{2,2} = length(Right(:,1)); B{2,3} = Right;
B{3,1} = 'Top';     B{3,2} = length(Top(:,1)); B{3,3} = flip(Top);
B{4,1} = 'Left';    B{4,2} = length(Left(:,1)); B{4,3} = flip(Left);

% Define the boundaries edges
bedges = cell(size(B,1),1);
for j = 1:size(B,1)    
    bedges{j} = bound(E,V,B{j,3});
end

% Get interior edge data
sizeS = max(max(E));
S = spalloc(sizeS,sizeS,6*length(E(:,1))); % 6 per cell 

% Create output matrix 
C = zeros(6*length(E(:,1)),4);
c = 1;

% Loop over all elements
for i = (1:length(E(:,1))) %i = triange number   
    for j = (1:3) 
        e = E(i,:); 
        e(j) = []; 
        e = sort(e); 
        
        if full(S(e(1),e(2))) ~= 0          
            t0 = full(S(e(1),e(2)));
            e0 = find(E(t0,:) ~= e(1) & E(t0,:) ~= e(2)); 
            if t0 < i %correctly order t1 and t2
                t1 = t0;
                e1 = e0;
                t2 = i;
                e2 = j;
            else
                t1 = i;
                e1 = j;
                t2 = t0;
                e2 = e0;
            end            
            C(c,:) = [t1 e1 t2 e2]; 
            c = c+1;            
        else
            S = S + sparse(e(1),e(2),i,sizeS,sizeS,6*length(E(:,1))); 
            %S(e(1),e(2)) = i            
        end
    end
end

C = C(1:c-1,:);
C = sortrows(C,[1 3]); 

%inedges = inedgedat(E,V,C);
% Now define the interior edge data for the mesh
inedges = zeros(size(C,1),8);

for i = 1:size(C,1)
    EL = C(i,1); 
    ER = C(i,3);
    tri = E(EL,:); %gives three node indicies of EL
    
    %conditions ensure that nA is always on top, and therefor normal points
    %towards R
    if C(i,2) == 1
        nA = tri(3); 
        nB = tri(2);
    elseif C(i,2) == 2
        nA = tri(1); 
        nB = tri(3);
    elseif C(i,2) == 3
        nA = tri(2); 
        nB = tri(1);
    end
    
    VA = V(nA,:); 
    VB = V(nB,:); %coordinates for both nodes  

    dl = ((VA(1)-VB(1))^2 + (VA(2)-VB(2))^2)^0.5;

    nx = (VA(2)-VB(2))/dl; 
    ny = (VB(1)-VA(1))/dl;
    
    inedges(i,1) = nA; %   nA is the index in V of node A
    inedges(i,2) = nB; %   nB is the index in V of node B
    inedges(i,3) = nx; %   nx and ny are the components of the unit normal vector respectively
    inedges(i,4) = ny; 
    inedges(i,5) = dl; %   dl is the length of the edge
    inedges(i,6) = EL; %   EL is the index of the left element in E
    inedges(i,7) = ER; %   ER is the index of the right element in E
end

end

function bedges = bound(E,V,B)

bedges = zeros(size(B,1),7); 

    for i = 1:size(B,1)
    
        nA = B(i,2); 
        nB = B(i,1);
        
        Eb = find(and(any(E == nA, 2), any(E == nB, 2)));   
        
        VA = V(nA,:); 
        VB = V(nB,:); 
        
        dl = ((VA(1)-VB(1))^2 + (VA(2)-VB(2))^2)^0.5;
        
        nx = (VA(2)-VB(2))/dl; 
        ny = (VB(1)-VA(1))/dl;
        
        bedges(i,1) = nA;
        bedges(i,2) = nB;
        bedges(i,3) = nx;
        bedges(i,4) = ny;
        bedges(i,5) = dl;
        bedges(i,6) = Eb;
      
    end
end

