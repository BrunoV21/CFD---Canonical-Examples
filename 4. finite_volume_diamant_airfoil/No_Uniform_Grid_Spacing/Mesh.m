function [x_export, y_export, nx, ny, dX_WW, dX_W, dX_E, dY_N, dY_S, indexfoil, cellsbefore, cellcount, x_airfoil] = Mesh
clc
%% Input for non uniform Mesh Generation

log_factor = 0.3; %input
DY0 = 0.005; %input
ymax = 4; %input, but it is not the final value
DY_max = 0.7; %input

x_airfoil = 2; %input
wake = 3; %input
pre_airfoil = 1;%input
xmax = 8;

Nx2 = 50; %Number of cells on the airofoil lenght

dX1_0 = 0.25; %Initial cell size before airfoil
GR1 = 0.9; %Growth rate

dX3_0 = 0.2; %Initial cell size after airfoil and wake
GR3 = 1.1;  %Growth rate

%% Function

function n = n(S, dX, GR)
    syms x 
    n = double(vpasolve( S == (dX*(1-GR^x))/(1-GR), x ));
end

%% Meshing

DY = DY0;
dY = [DY0];
Ny = 0; %Number of cells
y = 0; %aux 
while y <= ymax
   Ny = Ny + 1;
   y = y + DY;
   dY(end+1) = DY;
   if log(Ny + log_factor)*DY + DY0 < DY_max
      DY = log(Ny + log_factor)*DY + 1.5*DY0;
   end
      
end
Ny = Ny+1; %Number os rows in the mesh

%%%%%%%%%%%%%%%%%%5 x parameteress %%%%%%%%%%%%%%
x_airfoil_wake = x_airfoil + wake; 

x1 = x_airfoil - pre_airfoil;
x2 = x_airfoil_wake + pre_airfoil;
x3 = xmax;

n1 = floor(n(x1, dX1_0, GR1));
n3 = floor(n(xmax - x2, dX3_0, GR3));
dX1 = zeros(1, n1);
dX3 = zeros(1, n3);
dX1(1) = dX1_0;
dX3(1) = dX3_0;

for i = 2:n1
    dX1(i) = GR1*dX1(i-1); 
end


for i = 2:n3
   dX3(i) = GR3*dX3(i-1); 
end
sum(dX1);
sum(dX3);
dX2 = (xmax - sum(dX1) - sum(dX3))/Nx2 * ones(1, Nx2);
size(dX1);

dX = [dX1 dX2 dX3];

IDK = 10;

%Auxiliar
x = linspace(0.5*dX1_0, xmax-0.5*dX1_0, IDK);
%Horizontal lines
for i = 1:length(dY)
    plot(x, sum(dY(1:i))+zeros(IDK), 'black')
    hold on
end

ymax = sum(dY);
y = linspace(0.5*DY0, ymax, IDK);
%Vertical lines
for i = 1:length(dX)%length(dX)
    plot(sum(dX(1:i)) + zeros(IDK), y, 'black')
end

x_export = zeros(max(size(dX)),1);
x_export(1) = 0.5*dX1_0;
for i = 2:1: max(size(dX))
    x_export(i) = x_export(i-1) + dX(i);
end

y_export = zeros(max(size(dY)),1);
y_export(1) = 0.5*DY0;
for i = 2:1: max(size(dY))
    y_export(i) = y_export(i-1) + dY(i);
end

nx = max(size(x_export));
ny = max(size(y_export));

cellcount = nx*ny;

X = zeros(cellcount,1);

for j = 1:1:ny
    X((j-1)*nx+1:(j)*nx) = x_export(:);    
end

Y = zeros(cellcount,1);
y_export = flip(y_export);
for j = 1:1:nx
    Y((j-1)*ny+1:(j)*ny) = y_export(:);    
end

indexfoil = zeros(Nx2,1);
cellsbefore = max(size(dX1));
for i = 1:1:Nx2
    indexfoil(i) = cellsbefore +i;
end

%% rearranging the dx and dy vectors

% dX = 0.5 ;
% dY = dX;
% 
% xmax = 4;
% ymax = 0.75*xmax;
% 
% nx = xmax/dX;
% ny = ymax/dY;
% cellcount = nx*ny;


dX_WW = zeros(cellcount,1);
dX_W = zeros(cellcount,1);
dX_E = zeros(cellcount,1);

dY_N = zeros(cellcount,1);
dY_S = zeros(cellcount,1);
% 
% %% Interior Points
% BC_index = sort([1:1:nx, nx+1:nx:cellcount-nx, 2*nx:nx:cellcount-nx, cellcount-nx+1:1:cellcount]);

BC_index_x = sort([1, nx+1:nx:cellcount-nx, 2*nx:nx:cellcount-nx, cellcount]);
BC_index_xWB = sort([1, nx+1:nx:cellcount-nx]);
BC_index_xEB = sort([2*nx:nx:cellcount-nx, cellcount]); %East Boundary
BC_index_xW = sort([1, 2, nx+1:nx:cellcount-nx, nx+2:nx:cellcount-nx]);

BC_index_yN = 1:1:nx;
BC_index_yS = cellcount-nx+1:1:cellcount;

interior_points_X = setdiff(1:1:cellcount, BC_index_x);
interior_points_XW = setdiff(1:1:cellcount,BC_index_xW);

interior_points_YN = setdiff(1:1:cellcount, BC_index_yN);
interior_points_YS = setdiff(1:1:cellcount, BC_index_yS);

dX_WW(interior_points_XW) = X(interior_points_XW-1) - X(interior_points_XW-2);
dX_W(interior_points_X) = X(interior_points_X) - X(interior_points_X-1);
dX_E(interior_points_X) = X(interior_points_X+1) - X(interior_points_X);

dY_N(interior_points_YN) = Y(interior_points_YN-1) - Y(interior_points_YN);
dY_S(interior_points_YS) = Y(interior_points_YS) - Y(interior_points_YS+1);

dX_WW(BC_index_xW) = 0.5*dX1_0;
dX_W(BC_index_xWB) = 0.5*dX1_0;
dX_E(BC_index_xEB) = 0.5*dX3(end);

dY_N(1:1:nx) = 0.5*dY(end);
dY_S(cellcount-nx+1:1:cellcount) = 0.5*DY0;

%% Testting

% dX = 0.25;
% dY = dX;
% 
% xmax = 8;
% ymax = 0.75*xmax;
% 
% dX = 0.5 ;
% dY = dX;
% 
% xmax = 4;
% ymax = 0.75*xmax;
% 
% nx = xmax/dX;
% ny = ymax/dY;
% cellcount = nx*ny;
% 
% x_airfoil = 1;
% dY_N = zeros(cellcount, 1);
% dY_N(1:1:cellcount) = 0.5;
% dY_S = dY_N;
% dX_W =  dY_N;
% dX_WW =  dY_N;
% dX_E = dY_N;
% 
% dY_N(1:1:nx) = 0.25;
% dY_S(BC_index_yS) = 0.25;
% 
% dX_W(BC_index_xWB) = 0.25;
% dX_WW(BC_index_xWB) = 0;
% dX_WW(BC_index_xWB+1) = 0.25;
% 
% cellsbefore = round(x_airfoil/dX) ; %forcing the disance from w to the foil as 1m, same as the chord
% indexfoil = zeros(cellsbefore,1);
% for i = 1:1:cellsbefore
%     indexfoil(i) = cellcount-nx+i+cellsbefore;
% end
% 
% buildx = 1:1:nx;
% buildy = 1:1:ny;
% 
% %centroids :coo
% x_export = dX*0.5:buildx*dX:xmax;
% y_export = dY*0.5:buildy*dY:ymax;


end



