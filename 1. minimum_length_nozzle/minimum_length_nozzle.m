clc
close all
clear
warning('off')

fprintf('\n \n ============================ MINIMUM LENGTH NOZZLE ============================ \n\n\n ')

    while 1
            read = input('Choose Mach [1.3, 10]:\n\n ','s');
            if isnan(str2double(read)) ~= 1
                Ma_e = str2double(read);
                if Ma_e >= 1.3 && Ma_e <= 10 
                    break;
                end
            end
            fprintf( '\nWarning! Choose Mach contained in [1.3,10].' );
    end
    
    while 1
            readN = input('\nChoose the number of characteristics:\n\n ','s');
            if isnan(str2double(readN)) ~= 1
                N = str2double(readN);
                if N > 1 && (N-floor(N) == 0)
                    break;
                end
            end
            fprintf( '\n Warning! Choose a valid number of characteristics\n' );
    end
    
    while 1
      read2 = input('\nEnter 1 to continue with T_0 = 400, P_0 = 500 or 2 to change it:\n\n ');
           if read2 == 1 
                T_0 = 400;
                P_0 = 500;
                break;
           end
           if read2 == 2 
               while 1
                T_0 = input('\nT_0 = \n\n', 's');
                if isnan(str2double(T_0)) ~= 1
                T_0 = str2double(T_0);
                if T_0 > 0 
                    break;
                end
                
                end
                fprintf( '\nWarning! Enter a valid temperature' );
               end
               
                
               while 1 
                P_0 = input('P_0 = \n\n', 's');
                if isnan(str2double(P_0)) ~= 1
                P_0 = str2double(P_0);
                if P_0 > 0 
                    break;
                end
                
                end
                fprintf( '\nWarning! Enter a valid pressure' );
               end
               break
           end
    end

    fprintf('\n \n ============================ SOLVING ============================ \n\n ')

% Ma_e = 2.4

% N = 7; %number of divisions / iterations / characteristic

gamma = 1.4;


R = 287;

tic

Rho_0 = P_0/(T_0*R);

nu_e = (sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(Ma_e^2-1)))-atan(sqrt(Ma_e^2-1)))*180/pi; %Prandt Meyer Equation
theta_max = nu_e*0.5;


theta_0 = (theta_max-floor(theta_max));
dtheta = (theta_max - theta_0)/(N-1);


index = zeros (N, 1);
index(1) = N;
counts = index;
for i = 2:1:N
    counts(i) = N+1-i;
    index(i) = index(i-1)+N+1-i;
end

k_plus = zeros(index(end),1);
k_minus = zeros(index(end),1);
nu = zeros(index(end),1);
theta = zeros(index(end),1);
mach = zeros(index(end),1);
mu = zeros(index(end),1);

%1st section theta_0 -> theta_max // theta = nu
nu(1:(index(1))) = theta_0+dtheta*(0:1:(counts(1)-1));
theta = nu;
k_minus = 2*theta; % will be needed for the next steps


%2nd section theta from 0 to dtheta*index
[k_minus, k_plus] = N_section(k_plus, k_minus, index, counts);
theta = 0.5*(k_minus+k_plus);
nu = 0.5*(k_minus-k_plus);

for i=1:1:index(end)
[mach(i),nu(i),mu(i)] = flowprandtlmeyer(gamma,nu(i),'nu');
end

%% Plots

X = zeros(index(end),1);
Y = zeros(index(end),1);

% 1st line
X(1) = -cotd(theta(1)-mu(1));
Y(1) = 0;

%fist section
for i = 2:1:index(1)
    [X(i), Y(i)] = XandY(0, 1,X(i-1),Y(i-1),(theta(i)-mu(i)),0.5*(theta(i)+theta(i-1)+mu(i)+mu(i-1)));
end

for i = 2:1:N
    for j = 1:1:counts(i)
       point = index(i-1)+j;
       point_plus = index(i-1)+j-1;
       point_minus = index(i-1)-(counts(i)-j);
       S_plus = 0.5*(theta(point)+mu(point)+theta(point_plus)+mu(point_plus));
       S_minus = 0.5*(theta(point)-mu(point)+theta(point_minus)-mu(point_minus));
       if j == 1
           S_plus = 0;
           if i == 2
           point_plus = 1;
           else
           point_plus = index(i-2)+1;
           end
       end
       [X(point), Y(point)] = XandY(X(point_minus), Y(point_minus), X(point_plus), Y(point_plus), S_minus, S_plus);
     end  
end

X_nozzle = zeros(N,1);
Y_nozzle = zeros(N, 1);

for i = 1:1:index(1)
    if i == 1
        S_plus = 0.5*(theta(index(i))+mu(index(i))+theta(index(i)-1)+mu(index(i)-1));
        [X_nozzle(i), Y_nozzle(i)] = XandY(0, 1, X(index(i)), Y(index(i)), 0.5*(theta_max+theta(index(i))), S_plus);
    else
        S_plus = theta(index(i))+mu(index(i));
        [X_nozzle(i), Y_nozzle(i)] = XandY(X_nozzle(i-1), Y_nozzle(i-1), X(index(i)), Y(index(i)), 0.5*(theta(index(i-1))+theta(index(i))), S_plus);
    end
    
end    

Mach_centerline = zeros(index(1), 1);
Mach_nozzle = zeros(index(1), 1);
X_centerline = zeros(index(1), 1);
 for i = 1:1:index(1)
    if i == 1
        Mach_centerline(i) =  mach(1);
         X_centerline(i) = X(1);
    else
        Mach_centerline(i) =  mach(index(i-1)+1);
        X_centerline(i) = mach(index(i-1)+1);
    end
    Mach_nozzle(i) = mach(index(i));
 end
 

P_centerline = P_0./((1+(gamma-1)*0.5*Mach_centerline.^2).^(gamma/(gamma-1)));
P_nozzle = P_0./((1+(gamma-1)*0.5*Mach_nozzle.^2).^(gamma/(gamma-1)));

T_centerline = T_0./(1+0.5*(gamma-1)*Mach_centerline.^2);  
T_nozzle = T_0./(1+0.5*(gamma-1)*Mach_nozzle.^2);  

Rho_centerline = Rho_0./(1+0.5*(gamma-1)*Mach_centerline.^2).^(1/(gamma-1));
Rho_nozzle = Rho_0./(1+0.5*(gamma-1)*Mach_nozzle.^2).^(1/(gamma-1));

 toc
 fprintf('\n \n ============================ PLOTTING ============================ \n\n ')

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

tiledlayout(2,2)
title(sprintf('2D Minimum Lenght Nozzle'))
nexttile           

hold on
% axis equal
scatter([0, X_nozzle'], [1, Mach_nozzle'], 'red')
scatter([0, X_centerline'], [1, Mach_centerline'], 'blue')
plot([0 , X_nozzle(end)], [1, 1], 'black')
ylabel('Mach Number')
xlabel('lenght x/y0')
legend('Nozzle Contour','Centerline')
xlim([0 X_nozzle(end)]);
ylim([1, Mach_nozzle(end)]);

nexttile  

hold on
% axis equal
scatter([0, X_nozzle'], [P_0, P_nozzle'], 'red')
scatter([0, X_centerline'], [P_0, P_centerline'], 'blue')
plot([0 , X_nozzle(end)], [P_0, P_0], 'black')
ylabel('Pressure')
xlabel('lenght x/y0')
legend('Nozzle Contour','Centerline', 'P_0')
xlim([0 X_nozzle(end)]);
ylim([P_nozzle(end) P_0]);
nexttile  

hold on
% axis equal
scatter([0, X_nozzle'], [T_0, T_nozzle'], 'red')
scatter([0, X_centerline'], [T_0, T_centerline'], 'blue')
plot([0 , X_nozzle(end)], [T_0, T_0], 'black')
ylabel('Temperature')
xlabel('Lenght x/y0')
legend('Nozzle Contour','Centerline', 'T_0')
xlim([0 X_nozzle(end)]);
ylim([T_nozzle(end) T_0]);

nexttile  

hold on
% axis equal
scatter([0, X_nozzle'], [Rho_0, Rho_nozzle'], 'red')
scatter([0, X_centerline'], [Rho_0, Rho_centerline'], 'blue')
plot([0 , X_nozzle(end)], [Rho_0, Rho_0], 'black')
ylabel('Density')
xlabel('Lenght x/y0')
legend('Nozzle Contour','Centerline', 'Rho_0')
xlim([0 X_nozzle(end)]);
ylim([Rho_nozzle(end) Rho_0]);

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

hold on
axis equal
scatter(X, Y, '.', 'black')
scatter([0, X_nozzle'], [1, Y_nozzle'], 'red')

for i=1:1:index(1)
    if i == 1 
        %kPlus
        plot([X(1:1:index(i))', X_nozzle(i)'],[Y(1:1:index(i))', Y_nozzle(i)'], 'color', 'black') 
        %kMinus
        for j = 1:1:counts(i)
            plot([0, X(j)'],[1, Y(j)'], 'color', 'black')
            plot([0, X_nozzle(i)'], [1, Y_nozzle(i)'], 'color', 'blue')
        end
    else
        %kPlus
        plot([X(index(i-1)+1:1:index(i))', X_nozzle(i)'],[Y(index(i-1)+1:1:index(i))', Y_nozzle(i)'], 'color', 'black')
        %KMinus
        for j = 1:1:counts(i)
            point = index(i-1)+j;
            point_minus = index(i-1)-(counts(i)-j);
            plot([X(point_minus), X(point)],[Y(point_minus), Y(point)], 'color', 'black')
        end
        plot([X_nozzle(i-1)', X_nozzle(i)'], [Y_nozzle(i-1)', Y_nozzle(i)'], 'color', 'blue')
    end
end
xlim([0 X_nozzle(end)]);
ylim([0, Y_nozzle(end)]);
title(sprintf('2D Minimum Lenght Nozzle'))
         xlabel('Lenght x/y0')
         ylabel('Height y/y0')
         
 Ar = sqrt(1/(Ma_e)^2*(2/(gamma+1)*(1+0.5*(gamma-1)*Ma_e^2))^((gamma+1)/(gamma-1)));  
           
 fprintf('Area Ratio (relative to y0) = %.5f \n\n', Y_nozzle(end)) 
 fprintf('Length (relative to y0) = %.5f \n\n', X_nozzle(end)) 
 fprintf('Area Ratio from Isentropic Flow Tables(relative to y0) = %.5f \n\n', Ar)
 fprintf('Error relative to Isentropic Flow Area = %.2f %% \n\n', abs(Ar-Y_nozzle(end))/Ar*100)



%%
function [k_minus, k_plus] = N_section(k_plus, k_minus, index, counts)
    for i = 2:1: size(index)
        k_minus(index(i-1)+1:1:index(i)) = k_minus(i:1:counts(1));
        k_plus(index(i-1)+1:1:index(i)) = -k_minus(i);
    end
end

%%
function [X3,Y3] = XandY(X1,Y1,X2,Y2,Sminus1,Splus2)
% Point of intersection of two lines (given origin points and slopes)

X3 = (X1*tand(Sminus1)-X2*tand(Splus2)+(Y2-Y1))/(tand(Sminus1)-tand(Splus2));
       
Y3 = (tand(Sminus1)*tand(Splus2)*(X1-X2)+Y2*tand(Sminus1)-Y1*tand(Splus2))/...
                    (tand(Sminus1) - tand(Splus2));            
end
