close all
clc
clear

%% Physical Properties
gama = 1.4;
T0 = 288;
R = 287;


% Parameters Set up
% dX = 0.5 ;
% dY = dX;
% 
% xmax = 4;
% ymax = 0.75*xmax;
% 
% Mach = 1.8;
% 
% x_airfoil = 1;
%
% alpha = 0.65;
% 
% it_number = 30;
% 

%% Meshing
[x_export, y, nx, ny, dX_WW, dX_W, dX_E, dY_N, dY_S, indexfoil, cellsbefore, cellcount, x_airfoil] = Mesh;

%% Test Cases
% dX = 0.25;
% dY = dX;
% 
% xmax = 8;
% ymax = 0.75*xmax;
% 
% x_airfoil = 1;

Mach = 1.8;

alpha = 1.35;

it_number = 100;

%% Calculate Velocity

a_mach = sqrt(gama*R*T0);

U_inf = Mach*a_mach;

%% Shock Wave angles
beta_SA = asin(1/Mach); %Beta shock angle
beta_SA_degree = beta_SA * 180/pi;

delta_rad = atan(0.05); %Beta angle (wedge deflection)
delta_degree = 0.05 * 180/pi;

syms x
beta_OSA = vpasolve( 2/tan(x) * (Mach^2 *sin(x)^2 - 1)/((Mach^2 * (1.4 + cos(2*x))) + 2)  == 0.05, x, 0.6);  %Oblique Shock wave angle
beta_OSA_degree = beta_OSA * 180/3.1415926535;


%% Grid Initialization
% 
% nx = xmax/dX;
% ny = ymax/dY;
% cellcount = nx*ny;
% 
% %works as long as NX > NY
% buildx = 1:1:nx;
% buildy = 1:1:ny;
% 
% %centroids :coo
% x = dX*0.5:buildx*dX:xmax;
% y = dY*0.5:buildy*dY:ymax;


%% Airfoil Setup
% cellsbefore = round(x_airfoil/dX) ; %forcing the disance from w to the foil as 1m, same as the chord
% indexfoil = zeros(cellsbefore,1);

% for i = 1:1:cellsbefore
%     indexfoil(i) = cellcount-nx+i+cellsbefore;
% end


index_cells_positive_slope = indexfoil(1 : round(cellsbefore/2));
index_cells_negative_slope = indexfoil(round(cellsbefore/2)+1:end);

%% Solving
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
for line = index_cells_positive_slope(1):index_cells_positive_slope(end)
b(line) = 0.05*U_inf*(dX_W(line)+dX_WW(line))*dX_W(line)*dX_WW(line)*dY_N(line);
end
for line = index_cells_negative_slope(1):index_cells_negative_slope(end)
b(line) = -0.05*U_inf*(dX_W(line)+dX_WW(line))*dX_W(line)*dX_WW(line)*dY_N(line);
end

%% Building A Matrix
%% Boundary Condtions

%% Interior Points
for line = nx+3:1:cellcount-nx
    a(line, line) = ((1-Mach^2)*(dY_N(line)+dY_S(line))*dY_N(line)*dY_S(line)*dX_WW(line)-(dX_W(line)+dX_WW(line))*dX_W(line)*dX_WW(line)*(dY_N(line)+dY_S(line))); %P
    
    a(line, line-1) = -(1-Mach^2)*(dY_N(line)+dY_S(line))*dY_N(line)*dY_S(line)*(dX_WW(line)+dX_W(line)); %West
    a(line, line-2) = (1-Mach^2)*(dY_N(line)+dY_S(line))*dY_N(line)*dY_S(line)*(dX_W(line)); %West West
    
    a(line, line-nx) = (dX_W(line)+dX_WW(line))*dX_W(line)*dX_WW(line)*dY_S(line); % North
    a(line, line+nx) = (dX_W(line)+dX_WW(line))*dX_W(line)*dX_WW(line)*dY_N(line);   %South
end

a = forceBC_a(a, nx, cellcount, 0.5, 0.5, Mach)
% a(end,end)
% a(end-1,end-1)
% a(end-5,end-5)
% a(cellcount-1, cellcount-1)



%     format long

% a -

% [a_check, phi_check] = Desespero;
% 
% diagonala = zeros(cellcount,1);
% diagonalcheck = diagonala;
% for ii = 1:1:cellcount
%     diagonala(ii) = a(ii,ii);
%     diagonalcheck(ii) = a_check(ii,ii);
% end
% 2*diagonala- diagonalcheck %E
% a
% a_check
% a-a_check
% interior_points check P WW 
% North Inlet 1 e 2 e South
% South Inlet 1 e 2 e South
% phi_new-phi_check
%% Gauss Seidel and Jacobi Method
[phi_new, err] = Gauss_Siedel(a,b,phi, it_number, alpha, cellcount);
% [phi_new, err] = Jacobi(a,b,phi, it_number, alpha, cellcount);
% phi_new
% nx
% phi_new
% phi_check
%% Calculate U
%for interior points
for i = 2:1:cellcount-1
    u(i) = U_inf + (phi_new(i+1)-phi_new(i-1))/(dX_E(i)+dX_W(i));    
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
    v(i) = (phi_new(i+nx)-phi_new(i-nx))/(dY_N(i)+dY_S(i));
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


x_SW = linspace(x_airfoil, x(end)*3/4, 50);
b = @(beta) - tan(beta)*x_airfoil;
y_SW = @(x, beta, b) tan(beta) * x_SW + b + 0.5/2;

phi_plot = makesquare(phi_new, nx, ny, cellcount);
u_plot = makesquare(u, nx, ny, cellcount);
v_plot = makesquare(v, nx, ny, cellcount);

figure
hold on
title('Potential Field')
contourf(x_export,y,phi_plot, 500, 'LineColor','none')

figure
hold on
title('U - Velocity Field (Mach)')
contourf(x_export,y,u_plot/a_mach, 500, 'LineColor','none')
colormap(jet)
colorbar
% plot(x_SW, y_SW(x_SW, beta_SA, b(beta_SA)),'black')
% plot(x_SW, y_SW(x_SW, beta_OSA, b(beta_OSA)),'red')

figure
hold on
title('V - Velocity Field (m/s)')
contourf(x_export,y,v_plot, 500, 'LineColor','none')
colormap(jet)
colorbar

figure
semilogy(1:1:it_number, err)
         legend('Stream Function')
         xlabel('Iterations'); ylabel('Diference between iterations'); axis('square'); xlim([1 it_number]); grid on

fprintf("Shock Angle: %f\nOblique Shock Angle: %f\n", beta_SA_degree, beta_OSA_degree);



%% Functions
%% Force Boundary Conditions
function a = forceBC_a(a, nx, cellcount, dX, dY, Mach)
% first cell
a(1, :) = 0;
a(1,1) = 2*(1-Mach^2)*dY^2 - dX^2; %P
a(1,1+nx) = dX^2; %S

%2nd Corner
a(2, :) = 0;
a(2,2) = (1-Mach^2)*dY^2 - dX^2; %P
a(2,2+nx) = dX^2; %S
a(2,1) = -(1-Mach^2)*dY^2; %W


%3rd Corner
a(cellcount-nx+1, :) = 0;
a(cellcount-nx+1,cellcount-nx+1) = 2*(1-Mach^2)*dY^2 - dX^2; %P
a(cellcount-nx+1,cellcount-nx+2) = 0; %E
a(cellcount-nx+1,cellcount-nx+1-nx) = dX^2; %N


%4th Corner
a(cellcount-nx+2, :) = 0;
a(cellcount-nx+2,cellcount-nx+2) = (1-Mach^2)*dY^2 - dX^2; %P
a(cellcount-nx+2,cellcount-nx+1) =  -(1-Mach^2)*dY^2; %W
a(cellcount-nx+2,cellcount-nx+2-nx) = dX^2; %N


for cell = 3:1:nx %North BC - dPhi/dy = 0, 1st nx cells minus the first one
    a(cell, :) = 0;
    a(cell, cell) = (1-Mach^2)*dY^2 - dX^2; %P
    a(cell, cell + nx) = dX^2;  %S
    a(cell, cell-1) = -2*(1-Mach^2)*dY^2; %W
    a(cell, cell-2) = (1-Mach^2)*dY^2; %WW
end
    
for cell = cellcount-nx+3:1:cellcount %South BC - dPhi/dy = 0, last nx cells
    a(cell, :) = 0;
    a(cell, cell) = (1-Mach^2)*dY^2 - dX^2; %P
    a(cell, cell - nx) = dX^2;  %N
    a(cell, cell-1) = -2*(1-Mach^2)*dY^2; %W
    a(cell, cell-2) = (1-Mach^2)*dY^2; %WW
end

for cell = nx+1:nx:cellcount-nx %Inlet 1 BC - dPhi/dx = 0; and phi = 0 (first column)
    a(cell, :) = 0;
    a(cell, cell) = 2*((1-Mach^2)*dY^2-dX^2); %P
    a(cell, cell - nx) = dX^2;  %N
    a(cell, cell + nx) = dX^2;  %S   
end

for cell = nx+2:nx:cellcount-nx %Inlet 2 BC - dPhi/dx = 0 (second column)
    a(cell, :) = 0;
    a(cell, cell) = 3*(1-Mach^2)*dY^2 - dX^2; %P REVER MATH
    a(cell, cell - nx) = dX^2;  %N
    a(cell, cell + nx) = a(cell, cell - nx);  %South
    a(cell, cell-1) = -(1-Mach^2)*dY^2; %West;       
end


end


%% Gauss Seidel
function [phi_new, err] = Gauss_Siedel(A,b,phi,itnumber, alpha, cellcount)
phi_new = phi;
it = 1;
err = zeros(itnumber,1);
while it <= itnumber
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
function [phi_new, err] = Jacobi(A,b,phi,itnumber, alpha, cellcount)
phi_new = phi;
it = 1;
err = zeros(itnumber,1);

n = cellcount;

while it <= itnumber
    
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



