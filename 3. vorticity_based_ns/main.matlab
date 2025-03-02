function VorticityBasedNS
clear; close all


%Mesh Ref4 hi = 1.25cm  
% nx = 480; ny = 80;    
% ny2 = 120;  nx2 =960;
% dt = 0.001;

%Mesh Ref3 hi = 2.5cm
nx = 240; ny = 40;
ny2 = 60; nx2 =  480;
dt=0.001; 

%Mesh Ref2 hi = 5cm
% nx = 120; ny = 20;
% ny2 = 30; nx2 =  240;
% dt=0.001; 

%Mesh Ref1 hi = 10cm
% nx = 60;  ny = 10;
% ny2 = 15; nx2 = 120;
% dt=0.001; 

% Main Inputs
u_inlet = 1;
%Physics Continua
rho=1;                  %Density
nu=0.01;                %kinematic viscosity
tt=80;                  %total time
                        %time step

% nx = 60; ny = 10; 
Lx = 6; %retangulo 1
Ly = 1; %retangulo 1
dx = (Lx)/(nx); dxi = 1/dx;
dy = (Ly)/(ny); dyi = 1/dy;
x = linspace(-Lx,0,nx+1); 
y = linspace(0,Ly,ny+1);
%  xm = linspace(dx*0.5,0.5*dx+Lx,nx+1);
%  ym = linspace(dy*0.5,0.5*dy+Ly,ny+1);
[X,Y ]= meshgrid(x,y);


%Retangulo 2
% ny2 = 15;
% nx2 = 120;
Lx2 = 12; %retangulo 2
Ly2 = 1.5; %retangulo 2
dy2 = (Ly2)/(ny2); dyi2 = 1/dy2;
dx2 = (Lx2)/(nx2); dxi2 = 1/dx2;
x2 = linspace(0,Lx2,nx2+1); 
y2 = linspace(Ly-Ly2,Ly,ny2+1);
[X2,Y2] = meshgrid(x2,y2);

%BCs v_left (Inlet) // v_right (Outlet)
%    u_top (Top)   //  u_bottom (Ground)
v_left = 0; v_right = 0;
u_top = 0; u_bottom = 0;
% u_inlet = 1;


%Step
Lstep = 0.5;
nStep = round (Lstep/dx2) + 1;
xstream = 0.3*u_inlet;

%Physics Continua
% rho=1;                  %Density
% nu=0.01;                %kinematic viscosity
% tt=20;                   %total time
% dt=0.001;               %time step
ts = round(tt/dt);      %time step counter
% iits = 200; %round(0.5*tt/dt); %inner iterations counter
k = ts /800;             %graphic inner refresh rate

% block 1 matrix init
imin=2; imax=nx+imin-1;
jmin=2; jmax=ny+jmin-1;

% block 2 matrix init
imin2=2; imax2=nx2+imin2-1;
jmin2=2; jmax2=ny2+jmin2-1;

%minus         %Center       %plus   
im2 = 1:imax2-2; i2 = 2:imax2-1; ip2 = 3:imax2; 
jm2 = 1:jmax2-2; j2 = 2:jmax2-1; jp2 = 3:jmax2;
u2 = zeros(nx2+1,ny2+1);
un2 = zeros(nx2+1,ny2+1);
v2 = zeros(nx2+1,ny2+1);
vn2 = zeros(nx2+1,ny2+1);
% U = zeros(ny+1,nx+1); %Velocity (u, v)
% P = zeros(ny+1,nx+1); %Pressure goes away when creating
% Pn = zeros(ny+1,nx+1);% vorticity trasnport equations
% R = zeros(ny+1,nx+1);
W2 = zeros(nx2+1,ny2+1); %Vorticity
Wn2 = zeros(nx2+1,ny2+1); %NewVorticity used to store / update w
Q2 = zeros(nx2+1,ny2+1); %stream functon 
Qn2 = zeros(nx2+1,ny2+1);
WSS = zeros(nx2+1,ny2+1);
%Upwind 
%2
W2U = zeros(nx2+1,ny2+1);
%minus - 2         -1                 %Center           %plus +1    +2    %Upwind
im2U2 = 1:imax2-4; im2U1 = 2:imax2-3; i2U = 3:imax2-2;  ip2U1 = 4:imax2-1; ip2U2 = 5:imax2;
jm2U2 = 1:jmax2-4; jm2U1 = 2:jmax2-3; j2U = 3:jmax2-2;  jp2U1 = 4:jmax2-1; jp2U2 = 5:jmax2;
u2UMax = zeros(nx2+1,ny2+1);
u2UMin = zeros(nx2+1,ny2+1);
v2UMax = zeros(nx2+1,ny2+1);
v2UMin = zeros(nx2+1,ny2+1);
%1
W1U = zeros(nx+1,ny+1);
%minus - 2         -1                 %Center           %plus +1    +2    %Upwind
im1U2 = 1:imax-4; im1U1 = 2:imax-3; i1U = 3:imax-2;  ip1U1 = 4:imax-1; ip1U2 = 5:imax;
jm1U2 = 1:jmax-4; jm1U1 = 2:jmax-3; j1U = 3:jmax-2;  jp1U1 = 4:jmax-1; jp1U2 = 5:jmax;
u1UMax = zeros(nx+1,ny+1);
u1UMin = zeros(nx+1,ny+1);
v1UMax = zeros(nx+1,ny+1);
v1UMin = zeros(nx+1,ny+1);
% 
% size(u1UMax)

%vectorização
%minus         %Center       %plus    
im = 1:imax-2; i = 2:imax-1; ip = 3:imax; 
jm = 1:jmax-2; j = 2:jmax-1; jp = 3:jmax;
u = zeros(nx+1,ny+1);
v = zeros(nx+1,ny+1);
u_parabolic = zeros(nx2+1,ny2+1);
u_parabolic2 = zeros(nx2+1,ny2+1);
W = zeros(nx+1,ny+1); %Vorticity
Wn = zeros(nx+1,ny+1); %NewVorticity para guardar / atualizar w
Q = zeros(nx+1,ny+1); %stream functon 

erroW = zeros(1,ts+1);
erroQ = zeros(1,ts+1);
erroU = zeros(1,ts+1);
erroQ = zeros(1,ts+1);

Re = 1.5*u_inlet*(Ly2-Lstep)/nu;
courant = 1*dt*dxi;
fprintf('\n\nCourant number is %d\nRecomended CFL <<< 1.\n',courant)

fprintf('%d steps de %d s = %d s\n\nRunning...\n== ', ts, dt, tt)

u_inlet2 = u_inlet*Ly/Ly2;  

u_parabolic(1,j) =  6.*u_inlet.*(j-1)*dy/Ly.*(1-(j-1).*dy/Ly);  %Fully Developped Flow, u_inlet works as velocdiade média
u_parabolic2(1,j2) =  6.*u_inlet2.*(j2-1)*dy2/Ly2.*(1-(j2-1).*dy2/Ly2); %Fully Developped Flow at step apenas usado no fim para os graficos

folder = './'
for it = 1:ts+1 
    %Force BCs
    u(1,j) = u_parabolic(1,j); 
    
    
    %Inlet
    Q(1,j) = (4*Q(2,j)-Q(3,j)+2*dx*v(1,j))/3;
    W(1,j) = 2*dxi^2*(Q(1,j)-Q(2,j))-2*v(1,j)*dxi-0.5*dyi*(u(1,jp)-u(1,jm));
    %     W(1,1:jmax) = 2*dxi^2*(Q(1,1:jmax)-Q(2,1:jmax));
    
    %Outlet 
    Q(imax,j) = (4*Q(imax-1,j)-Q(imax-2,j)-2*dx*v(imax,j))/3;
    W(imax,j) = 2*dxi^2*(Q(imax,j)-Q(imax-1,j))+2*dxi*v(imax,j)-0.5*dyi*(u(imax,jp)-u(imax,jm));
%     Q(imax,1:jmax) = 0;
%     W(imax,1:jmax) = 2*dxi^2*(Q(imax,1:jmax)-Q(imax-1,1:jmax));
        
    %Top
    Q(1:imax,jmax) = u_inlet*Ly; %Check this later 
    W(1:imax,jmax) = 2*dyi^2*(Q(1:imax,jmax)-Q(1:imax,jmax-1))-2*dyi*u_top;
    
    %Bottom
    Q(1:imax,1) = 0;
    W(1:imax,1) = 2*dyi^2*(Q(1:imax,1)-Q(1:imax,2))+2*dyi*u_bottom;
   

    
    %Vorticity Transport Eq
    Wn = W;
    
        %Upwind Crl1
%     arrayfun(@upwindMax, u2) arrayfun(@upwindMin, u2)
    u1UMax(i,j) = 0.5*(u(i,j)+abs(u(i,j)));
    u1UMin(i,j) = 0.5*(u(i,j)-abs(u(i,j)));
    v1UMax(i,j) = 0.5*(v(i,j)+abs(v(i,j)));
    v1UMin(i,j) = 0.5*(v(i,j)-abs(v(i,j)));
    
    W1U(i1U,j1U) = (...
          u1UMax(i1U,j1U).*(W(im1U2,j1U)-3*W(im1U1,j1U)+3*W(i1U,j1U)-W(ip1U1,j1U))*dxi/3+...
          u1UMin(i1U,j1U).*(W(im1U1,j1U)-3*W(i1U,j1U)+3*W(ip1U1,j1U)-W(ip1U2,j1U))*dxi/3+...
          v1UMax(i1U,j1U).*(W(i1U,jm1U2)-3*W(i1U,jm1U1)+3*W(i1U,j1U)-W(i1U,jp1U1))*dyi/3+...
          v1UMin(i1U,j1U).*(W(i1U,jm1U1)-3*W(i1U,j1U)+3*W(i1U,jp1U1)-W(i1U,jp1U2))*dyi/3);
    
    W(i,j) = Wn(i,j)+...
        dt*(-(Q(i,jp)-Q(i,jm))*0.5*dyi.*(Wn(ip,j)-Wn(im,j))*0.5*dxi+...
        (Q(ip,j)-Q(im,j))*0.5*dxi.*(Wn(i,jp)-Wn(i,jm))*0.5*dyi+...
        nu*dxi*dyi*(Wn(ip,j)+Wn(i,jp)-4*Wn(i,j)+Wn(im,j)+Wn(i,jm))-...
        min(1.2*dt*max(max(abs(u(i,j)))*dxi,max(abs(v(i,j)))*dyi),1).*W1U(i,j)); %doublecheck last term denominador
    %Get Stream Function Q
    Q(i,j) = 0.25*(W(i,j)*dx*dy+Q(ip,j)+Q(i,jp)+Q(i,jm)+Q(im,j));
    %Get u e v from Q
%     u(1,2:jmax-1) = u_parabolic(1,j);
%     u(imax,1:jmax) = u(imax-1,1:jmax);
    u(1:imax,jmax) = u_top;
    u(1:imax,1) = u_bottom;
    u(i,j) = 0.5*dyi*(Q(i,jp)-Q(i,jm));
    v(i,j) = -0.5*dxi*(Q(ip,j)-Q(im,j));
    u(imax,j) = u(imax-1,j);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2
    %Force BCs
    u2(1,nStep+1:jmax2-1) = u(imax,j);
    u2(imax2,j2) = u2(imax2-1,j2);
    
    %Inlet
    Q2(1,nStep+1:jmax2-1) = Q(imax, j);
    W2(1,nStep+1:jmax2-1) = W(imax, j);
    W2(1,1:nStep) = -2*Q2(2,1:nStep) /(dx2^2);  %step
    %     W(1,1:jmax) = 2*dxi^2*(Q(1,1:jmax)-Q(2,1:jmax));
    
    %Outlet 
    Q2(imax2,j2) = (4*Q2(imax2-1,j2)-Q2(imax2-2,j2)-2*dx2*v2(imax2,j2))/3;
    W2(imax2,j2) = 2*dxi2^2*(Q2(imax2,j2)-Q2(imax2-1,j2))+2*dxi2*v2(imax2,j2)-0.5*dyi2*(u2(imax2,jp2)-u2(imax2,jm2));
%     Q(imax,1:jmax) = 0;
%     W(imax,1:jmax) = 2*dxi^2*(Q(imax,1:jmax)-Q(imax-1,1:jmax));
        
    %Top
    Q2(1:imax2,jmax2) = u_inlet*Ly; %Check this later 
    W2(1:imax2,jmax2) = 2*dyi2^2*(Q2(1:imax2,jmax2)-Q2(1:imax2,jmax2-1))-2*dyi2*u_top;
    
    %Bottom
    Q2(1:imax2,1) = 0;
    W2(1:imax2,1) = 2*dyi2^2*(Q2(1:imax2,1)-Q2(1:imax2,2))+2*dyi2*u_bottom;
    
    
    %Vorticity transport Eq Block 2
    Wn2 = W2;
    
    %Upwind Crl2
%     min(1.2*dt*max(max(u(i,j),max(v(i,j))),1))
%     arrayfun(@upwindMax, u2) arrayfun(@upwindMin, u2)
    u2UMax(i2,j2) = 0.5*(u2(i2,j2)+abs(u2(i2,j2)));
    u2UMin(i2,j2) = 0.5*(u2(i2,j2)-abs(u2(i2,j2)));
    v2UMax(i2,j2) = 0.5*(v2(i2,j2)+abs(v2(i2,j2)));
    v2UMin(i2,j2) = 0.5*(v2(i2,j2)-abs(v2(i2,j2)));
    W2U(i2U,j2U) = ...
          u2UMax(i2U,j2U).*(W2(im2U2,j2U)-3*W2(im2U1,j2U)+3*W2(i2U,j2U)-W2(ip2U1,j2U))*dxi2/3+...
          u2UMin(i2U,j2U).*(W2(im2U1,j2U)-3*W2(i2U,j2U)+3*W2(ip2U1,j2U)-W2(ip2U2,j2U))*dxi2/3+...
          v2UMax(i2U,j2U).*(W2(i2U,jm2U2)-3*W2(i2U,jm2U1)+3*W2(i2U,j2U)-W2(i2U,jp2U1))*dyi2/3+...
          v2UMin(i2U,j2U).*(W2(i2U,jm2U1)-3*W2(i2U,j2U)+3*W2(i2U,jp2U1)-W2(i2U,jp2U2))*dyi2/3;
    
    
    
    W2(i2,j2) = Wn2(i2,j2)+...
        dt*(-(Q2(i2,jp2)-Q2(i2,jm2))*0.5*dyi2.*(Wn2(ip2,j2)-Wn2(im2,j2))*0.5*dxi2+...
        (Q2(ip2,j2)-Q2(im2,j2))*0.5*dxi2.*(Wn2(i2,jp2)-Wn2(i2,jm2))*0.5*dyi2+...
        nu*dxi2*dyi2*(Wn2(ip2,j2)+Wn2(i2,jp2)-4*Wn2(i2,j2)+Wn2(im2,j2)+Wn2(i2,jm2))-...
        min(1.2*dt*max(max(abs(u2(i2,j2)))*dxi2,max(abs(v2(i2,j2)))*dyi2),1).*W2U(i2,j2)); %doublecheck no denominador do ultimo termo
    %Get Stream Function Q
    Qn2 = Q2;
    Q2(i2,j2) = 0.25*(W2(i2,j2)*dx2*dy2+Q2(ip2,j2)+Q2(i2,jp2)+Q2(i2,jm2)+Q2(im2,j2));
    un2 = u2;
    vn2 = v2;
    u2(1:imax2,jmax2) = u_top;
    u2(1:imax2,1) = u_bottom;
    u2(i2,j2) = 0.5*dyi2*(Q2(i2,jp2)-Q2(i2,jm2));
    v2(i2,j2) = -0.5*dxi2*(Q2(ip2,j2)-Q2(im2,j2));
    u2(imax2,j2) = u2(imax2-1,j2);
    
    %Erro /Convergence
    erroW(1,it) = max(max(W2-Wn2));
    erroQ(1,it) = max(max(Q2-Qn2));
    erroU(1,it) = max(max(u2-un2)); 
    erroV(1,it) = max(max(v2-vn2));
        
       if (it/k) == floor(it/k) %Plot
         %Stream Function
         f1 = figure(1);
         f1.Position = [5 525 2095 250]; 
         %% 
%          clf
         hold on
         t = tiledlayout(1,3);
         starty = [-0.15, -0.2, -0.25, -0.3, -0.35, -0.4, -0.45, -0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
         startx = [xstream, xstream, xstream, xstream, xstream, xstream, xstream, xstream, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
         startxCh = -Lx*zeros(size(starty));
         ax1 = nexttile;
         h1 = streamline(x,y,-u',-v',startxCh,starty,[1, 1000]);
%          area(-Lx, Ly-Ly2)
         axis('equal',[-Lx -dx 0 Ly]); 
         set(h1,'color','k')
         patch([-Lx 0 0 -Lx],[0 0 Ly-Ly2 Ly-Ly2],'black')
         ax2 = nexttile([1 2]);
%          size(starty)
%          size(startx)
         %0.3*ones(size(starty));%max(x2)*rand(1000,1); 
%        starty = [0, 0.001, 0.01, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5]
         h = streamline(x2,y2,u2',v2',startx,starty,[1, 1000]);
         axis('equal',[0 Lx2 Ly-Ly2 Ly]); 
         set(h,'color','k')
         
         linkaxes([ax1,ax2],'y');
         title(t,sprintf('Stream Function\nRe = %d',Re))
         xlabel(t,'x-location')
         ylabel(t,'y-location')
         yticklabels(ax2,{})
         t.TileSpacing = 'none';
         hold off
         pngFileName = sprintf('plot_subject_%d.png', it);
         fullFileName = fullfile(folder, pngFileName);
         fullFileName = fullfile(folder, pngFileName);
         saveas(gcf, fullFileName);
%          exportgraphics(gcf, sprintf('Fig_%d.png', it));
         drawnow
         
                 
%          U Magnitude
%          f2 = figure(2);
%          f2.Position = [475 250 200 250]; 
%          clf
%          hold on
%          set(h1,'color','k')
%          patch([-Lx 0 0 -Lx],[0 0 Ly-Ly2 Ly-Ly2],'black')
%          contourf(X,Y,W',23,'LineColor','none');
%          quiver(x, y, u', v',.7,'k-')
%          title(sprintf('U - Magnitude\nRe = %d   t = %.3g s',Re,it*dt))
%          xlabel('x-location')
%          ylabel('y-location')
%          axis('equal',[-Lx 0 Ly-Ly2 Ly]); 
%          colormap jet
%          colorbar
%          hold off+
%          drawnow
%          Pressure Scalar + Velocity Vector
%          f3 = figure(3);
%          f3.Position = [1005 250 600 525]; 
%          clf
%          hold on
%          contourf(X2,Y2,W2',23,'LineColor','none');
%          quiver(x2, y2, u2', v2',.7,'k-')
%          title(sprintf('V-Magnitude\nRe = %d   t = %.3g s',Re,it*dt))
%          xlabel('x-location')
%          ylabel('y-location')
%          axis('equal',[0 Lx2 Ly-Ly2 Ly]);
%          colormap jet
%          colorbar
%          hold off
%          drawnow
%          fprintf('%d%% == ',it/ts*100)
       end
 end
%          Proof
%        Vorticity at bottom
         figure(4); plot(x2, W2(:,1));
         title('Vorticity at bottom');
         line([0 Lx2],[-u_parabolic2(1,2)*dyi2 -u_parabolic2(1,2)*dyi2],'Color','black','LineStyle','--')
         xlabel('U'); ylabel('y'); axis('square'); xlim([0 Lx2]); grid on
%        Vorticity at top 
         figure(5); plot(x2,W2(:,jmax2));
         title('Vorticity a top');
         line([0 Lx2],[u_parabolic2(1,jmax2-1)*dyi2 u_parabolic2(1,jmax2-1)*dyi2],'Color','black','LineStyle','--')
         xlabel('x'); ylabel('W'); axis('square'); xlim([0 Lx2]); grid on
% Velocity at Step   
         figure(6); 
         x1 = u_parabolic(1,1:jmax); y1 = (0:jmax-1)*dy/Ly;
         x2 = u2(1,nStep+1:jmax2-1); y2 =(j-1)*dy/Ly;
         plot(x1,y1, x2, y2, 'o')
      

         title('Velocity at Step');
         xlabel('U'); ylabel('y'); axis('square'); xlim([0 1.5*u_inlet]); grid on

% Velocity at Outlet   
         figure(7); 
         x1 = u_parabolic2(1,1:jmax2); y1 = -0.5:dy:1;
         x2 = u2(imax2,1:jmax2); y2 = -0.5:dy:1;
         plot(x1,y1, x2, y2, 'o')

         title('Velocity at Outlet');
         xlabel('U'); ylabel('y'); axis('square'); xlim([0 1.5*u_inlet]); grid on

% Erro 
         figure(8); 
         x = 1:ts+1; y1 = erroW(1,1:ts+1);
                     y2 = erroQ(1,1:ts+1);
                     y3 = erroU(1,1:ts+1);
                     y4 = erroV(1,1:ts+1);
%         Azul Vorticidade %Vermelhor Strean Function
         semilogy(x, y1, x, y2, x, y3, x, y4)
         legend('Vorticity','Stream Function', 'U-Velocity', 'V-Velocity')
         xlabel('Número de Iterações'); ylabel('Diferença entre iterações'); axis('square'); xlim([1 ts+1]); grid on
         
                  


% min(1.2*dt*max(max(abs(u2(i2,j2)))*dxi2,max(abs(v2(i2,j2)))*dyi2),1).*W2U(i2,j2)
%   u2UMax(i2,j2) = 0.5*(u2(i2,j2)+abs(u2(i2,j2)))
%     u2UMin(i2,j2) = 0.5*(u2(i2,j2)-abs(u2(i2,j2)))
%     v2UMax(i2,j2) = 0.5*(v2(i2,j2)+abs(v2(i2,j2)))
%     v2UMin(i2,j2) = 0.5*(v2(i2,j2)-abs(v2(i2,j2)))

WSS = -W2.*nu;             
for length = 1:imax2         
  if WSS(length,1) > 0
      if length*dx2 > 5*dx2
      fprintf('\n\nCalculations Completed\n\nReattachment length = %.2dm',length*dx2)
      break
      end
   end
 end      
fprintf('\nCalculations Completed\n\nEnjoy the the Coulourful Fluid Diagrams :)\n\n')

end