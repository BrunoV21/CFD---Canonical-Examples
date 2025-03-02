
close all
clc
clear

%% Physical Properties
gama = 1.4;
T0 = 288;
R = 287;
p_atm = 1.01325e5;
rho_0 = p_atm/(R*T0);

%% Parameters Set up
% dX = 0.5 ;
% dY = dX;
% 
% xmax = 4;
% ymax = 0.75*xmax;
% 
% x_airfoil = 1;
%
% Mach = 0.85;
% 
% alpha = .65;
% 
% it_number = 30;


%% Test Cases
dX = 0.25;
dY = dX;

xmax = 8;
ymax = 0.75*xmax;

x_airfoil = 2.5;

Mach = 0.95;

alpha = 1.3;

it_number = 2000;

%% Calculate Velocity

a_mach = sqrt(gama*R*T0);

U_inf = Mach*a_mach;

%% Grid Initialization

nx = xmax/dX;
ny = ymax/dY;
cellcount = nx*ny;

%works as long as NX > NY
buildx = 1:1:nx;
buildy = 1:1:ny;

%centroids :coo
x = dX*0.5:buildx*dX:xmax;
y = dY*0.5:buildy*dY:ymax;

%% Airfoil Setup
cellsbefore = round(x_airfoil/dX) ; %forcing the disance from w to the foil as 1m, same as the chord
indexfoil = zeros(cellsbefore,1);

for i = 1:1:cellsbefore
    indexfoil(i) = cellcount-nx+i+cellsbefore;
end


index_cells_positive_slope = indexfoil(1 : round(cellsbefore/2));
index_cells_negative_slope = indexfoil(round(cellsbefore/2)+1:end);

tic
%% Matrix Initialiation
a = zeros(cellcount, cellcount);
phi = zeros(cellcount, 1);
phi_new = phi;
T1 = a;
T2 = a;
%% Coefficients initialization
b = phi; %sources
u = phi;
v = phi;

%% Interior Points
BC_index = sort([1:1:nx, nx+1:nx:cellcount-nx, 2*nx:nx:cellcount-nx, cellcount-nx+1:1:cellcount]);
interior_points = setdiff(1:1:cellcount,BC_index);

%% Boundary Condtions
a = forceBC_a(a, nx, cellcount, dX, dY, Mach);

%% Gauss Seidel
[phi_new, err] = Solve(a,b,phi,it_number, alpha, dX, dY, Mach, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, gama ,nx, cellcount);


%% Calculate U
%for interior points
for i = 2:1:cellcount-1
    u(i) = U_inf + 0.5*(phi_new(i+1)-phi_new(i-1))/dX;    
end
for i = nx:nx:cellcount
    u(i) = u(i-1);
end
for i = 1:nx:cellcount
    u(i) = U_inf;
end
% u(indexfoil) = 0; no slip condition

%% Calculate V
for i = nx+1:1:cellcount-nx
    v(i) = 0.5*(phi_new(i+nx)-phi_new(i-nx))/dX;
end
for i = 1:1:nx
    v(i) = v(i+nx);
end
for i = (ny-1)*nx:1:cellcount
   v(i) = 0;
end
for i = index_cells_positive_slope
    v(i) = -0.05*U_inf;
end
for i = index_cells_negative_slope
    v(i) = 0.05*U_inf;
end
% v(indexfoil) = 0; impermeability condityion

toc
%% Prepare a variable to be plotted
phi_plot = makesquare(phi_new, nx, ny, cellcount);
u_plot = makesquare(u, nx, ny, cellcount);
v_plot = makesquare(v, nx, ny, cellcount);

figure
hold on
title('Potential Field')
contourf(x,y,phi_plot, 500, 'LineColor','none')

figure
hold on
title('U - Velocity Field (Mach)')
contourf(x,y,u_plot/a_mach, 500, 'LineColor','none')
colormap(jet)
colorbar

figure
hold on
title('V - Velocity Field (m/s)')
contourf(x,y,v_plot, 500, 'LineColor','none')
colormap(jet)
colorbar

rho = zeros(cellcount,1);
for i = 1:1:cellcount
    rho(i) = rho_0*(1+0.5*(gama-1)*(u(i)/a_mach)^2)^(-1/(gama-1));
end

rho_plot = makesquare(rho, nx, ny, cellcount);

figure
hold on
title('Rho - Density Field (kg/m^3)')
contourf(x,y,rho_plot, 500, 'LineColor','none')
colormap(jet)
colorbar

figure
semilogy(1:1:it_number, err)
         legend('Stream Function')
         xlabel('Iterations'); ylabel('Diference between iterations'); axis('square'); xlim([1 it_number]); grid on



%% Functions
%% Force Boundary Conditions
function a = forceBC_a(a, nx, cellcount, dX, dY, Mach)
% first cell
a(1, :) = 0;
a(1,1) = 3*dY^2*(1-Mach^2)+dX^2; %P
a(1,2) = - -(1-Mach^2)*dY^2; %E
a(1,1+nx) = -dX^2; %S

%2nd Corncer
a(nx, :) = 0;
a(nx,nx) = dY^2*(1-Mach^2)+dX^2; %P
a(nx,nx-1) = - -(1-Mach^2)*dY^2; %W
a(nx,nx+nx) = -dX^2; %S


%3rd Corner
a(cellcount-nx+1, :) = 0;
a(cellcount-nx+1,cellcount-nx+1) = 3*dY^2*(1-Mach^2)+dX^2; %P
a(cellcount-nx+1,cellcount-nx+2) = - -(1-Mach^2)*dY^2; %E
a(cellcount-nx+1,cellcount-nx+1-nx) = -dX^2; %N


%Last cell
a(cellcount, :) = 0;
a(cellcount,cellcount) = dY^2*(1-Mach^2)+dX^2; %P
a(cellcount,cellcount-1) = - -(1-Mach^2)*dY^2; %W
a(cellcount,cellcount-nx) = -dX^2; %N




for cell = 2:1:nx-1 %North BC - dPhi/dy = 0, 1st nx cells minus the first one
    a(cell, :) = 0;
    a(cell, cell) = 2*dY^2*(1-Mach^2)+dX^2; %P
    a(cell, cell + nx) = -dX^2;  %South
    a(cell, cell-1) = -(1-Mach^2)*dY^2; %West
    a(cell, cell+1) = a(cell, cell+1); %East
end
    
for cell = cellcount-nx+2:1:cellcount-1 %South BC - dPhi/dy = 0, last nx cells
    a(cell, :) = 0;
    a(cell, cell) = 2*dY^2*(1-Mach^2)+dX^2; %P
    a(cell, cell - nx) = -dX^2;  %North
    a(cell, cell-1) = -(1-Mach^2)*dY^2; %West
    a(cell, cell+1) = a(cell, cell-1); %East
end

for cell = nx+1:nx:cellcount-nx %West BC - dPhi/dx = 0; and phi = 0 
    a(cell, :) = 0;
    a(cell, cell) = 3*dY^2*(1-Mach^2)+2*dX^2; %P
    a(cell, cell - nx) = -dX^2;  %North
    a(cell, cell + nx) = a(cell, cell - nx);  %South
    a(cell, cell+1) = -(1-Mach^2)*dY^2; %East
end

for cell = 2*nx:nx:cellcount-nx %East BC - dPhi/dx = 0; 
    a(cell, :) = 0;
    a(cell, cell) = 2*dY^2*(1-Mach^2)+2*dX^2; %P
    a(cell, cell - nx) = -dX^2;  %North
    a(cell, cell + nx) = a(cell, cell - nx);  %South
    a(cell, cell-1) = -(1-Mach^2)*dY^2; %West;       
end


end

%% Calculate Z
function Z = calculateZ(Mach, U_inf, gama, dX, phi, interior_points)
Z(interior_points) = 1 - Mach^2 - (1+gama)*Mach^2/U_inf*(phi(interior_points-1) - phi(interior_points-2))/dX;
end

%% Build A and B
function [a, b] = buildA_B(a, b, dX, dY, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, nx, Z)
for line = interior_points
    
    a(line, line) = max(Z(line),0)/Z(line)*2*(dX^2+dY^2*Z(line)) + min(Z(line),0)/Z(line)*(dY^2*Z(line)-2*dX^2);
    a(line, line-nx) = min(Z(line),0)/Z(line)*dX^2+max(Z(line),0)/Z(line)*(-dX^2);
    a(line, line-1) = min(Z(line),0)*-2*dY^2+max(Z(line),0)*-dY^2;
    a(line, line-2) = min(Z(line),0)*dY^2;
    a(line, line+1) =  max(Z(line),0)*-dY^2;
    a(line, line+nx) = min(Z(line),0)/Z(line)*dX^2+max(Z(line),0)/Z(line)*(-dX^2);   
    
end
b(index_cells_positive_slope(1):index_cells_positive_slope(end)) = max(Z,0)/Z*-0.05*dX^2*dY*U_inf + min(Z,0)/Z*0.05*dX^2*dY*U_inf;
b(index_cells_negative_slope(1):index_cells_negative_slope(end)) = -(max(Z,0)/Z*-0.05*dX^2*dY*U_inf + min(Z,0)/Z*0.05*dX^2*dY*U_inf);
end

%% Gauss Seidel
function [phi_new, err] = Solve(A,b,phi,itnumber, alpha, dX, dY, Mach, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, gama ,nx, cellcount)
phi_new = phi;
it = 1;
err = zeros(itnumber,1);
Z = phi;
while it <= itnumber
    Z = calculateZ(Mach, U_inf, gama, dX, phi, interior_points);
    [A, b] = buildA_B(A, b, dX, dY, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, nx, Z);
    for i = 1:1:cellcount
        Sum_new = 0;
        Sum = 0;
        for j = 1:1:i-1
            Sum_new = Sum_new + A(i, j)*phi_new(j);
        end
        for j = i+1:1:cellcount
            Sum = Sum + A(i, j)*phi(j);
        end        
        phi_new(i) = (1-alpha)*phi(i)+alpha/(A(i,i))*(b(i)- Sum_new - Sum);
    end
    err(it) = sum(abs(phi_new - phi));                    % finding error
    it = it + 1;
    phi = phi_new;
end
end

%% Jacobi
function [phi_new, err] = Jacobi(A,b,phi,itnumber, alpha, dX, dY, Mach, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, gama ,nx, cellcount)
phi_new = phi;
it = 1;
err = zeros(itnumber,1);

n = cellcount;

while it <= itnumber
    Z = calculateZ(Mach, U_inf, gama, dX, phi, interior_points);
    [A, b] = buildA_B(A, b, dX, dY, interior_points, index_cells_positive_slope, index_cells_negative_slope, U_inf, nx, Z);
    for i=1:n
        sigma=0;
        for j=1:n  
            if j~=i
                sigma=sigma+A(i,j)*phi(j);
            end     
        end     
        
        phi_new(i) = (1/A(i,i))*(b(i)-sigma);
    end
    
    err(it) = sum(abs(phi_new - phi));                    % finding error
    it = it + 1;
    phi = phi_new;
end

end

%% Prepare a variable to be plotted
function phi_plot = makesquare(phi, nx, ny, cellcount)
phi_plot = zeros(ny, nx);
for i = ny:-1:1
    phi_plot(i,:) = phi((ny-i)*nx+1 : (ny-i+1)*nx);
end
end