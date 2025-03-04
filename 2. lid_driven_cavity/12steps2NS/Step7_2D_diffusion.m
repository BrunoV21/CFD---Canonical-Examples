%%%variable declarations
clear all
nx = 31;
ny = 31;
nt = 170;
dx = 2.0/(nx-1);
dy = 2.0/(ny-1);
sigma = .25;
nu=.05
dt = sigma*dx*dy/nu;

x = linspace(0,2,nx);
y = linspace(0,2,ny);

u = ones(ny,nx); %%create a 1xn vector of 1's
un=ones(ny,nx);

%%%Assign initial conditions

u(.5/dy:1/dy+1,.5/dx:1/dx+1)=2;


%%%Plot Initial Condition
 %%the figsize parameter can be used to produce different sized images
                  
[X, Y] = meshgrid(x,y);                            


for n=1:nt+1
    un=u;
    for i=2:(ny-1)
        for j=2:(nx-1)
        u(i,j)=un(i,j)+nu*dt/dx^2*(un(i+1,j)-2*un(i,j)+un(i-1,j))+nu*dt/dy^2*(un(i,j+1)-2*un(i,j)+un(i,j-1));
        end
    end
   surf(x,y,u)
   pause(0.1)
end

