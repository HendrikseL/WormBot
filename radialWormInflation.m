%{
Calculates the pressure after a set itme of inflatio for the end soft actuators.
%}

clear
close all
tic

%Fluid Properties
Ptank = 101325;%starts at atmo pressure 
R = 287; %J/kgK
Pcomp = 10000 + Ptank; %estimated compressor pressure 
Tstag = 298; %K, isothermal model Tstag will not change
rhoTank = Ptank/(Tstag*R);
gamma = 1.4;

%geometric properties
r_init = (75/2)*1e-3; %m
l_init = 0.075; %initial length
thickness = 5e-3; %m
E = 69e3; %young's modulus, Pa
Avalve = 3.1416e-06; %area of the valve, circle of radius 10mm

%stagnation density of fluid in the compressor
rhoComp= ((Pcomp/Ptank)^(1/gamma))*rhoTank;
%inital volume of fully inflated tank
[volume r] = updateVolume(Ptank, l_init,r_init,thickness, E); %m^3 

%initial mass
m = rhoTank * volume;

%under relaxation factor
a = 0.01;

%from worm_inflate, mass flow of the compressor
mdot = 0.00048188; 

deltat = 0.0001;
t = 0;
t_target = 0.8;
i = 1;
%main loop
while t < t_target
   
    %update stagnation density and total mass
    m(i+1) = m(i) + (mdot * deltat);
    rhoTank = m(i+1) / volume(i);
    
   
    %inner loop, checking for volume convergence
    flag = 1;
    j = 1;
    vol_inter = zeros(1000,1);
    vol_inter(1) = volume(i);
    while flag == 1 && j < 4000
        % update stagnation pressure
        Ptank(i+1) = Pcomp*(rhoTank/rhoComp)^1.4;
    
        %update volume
        [vol_update r_update]= updateVolume(Ptank(i+1), l_init,r_init,thickness, E);
        vol_inter(j+1) = (1-a)*vol_inter(j) + a*vol_update;

        rhoTank = m(i+1) / vol_inter(j+1);

        %convergence check
        if abs(vol_inter(j+1)-vol_inter(j)) < 1e-7
            flag = 0;
            volume(i+1) = vol_inter(j+1);
            r(i+1) = r_update;
        end

        %Convergence iteration
        j = j +1;
    end

    %Update mach number
    M(i+1) = real(sqrt( (2/(gamma-1)) * ( (Pcomp/Ptank(i+1))^((gamma-1)/gamma) -1 )));

    if M(i+1) > 1
        M(i+1) =1;
    end


    t(i+1) = t(i) + deltat;
    i = i + 1;
end

fprintf("Discharge Time: %fl seconds\n",t(end));
toc
volumeEnd = pi*r_init^2*l_init;
figure(1)
hold on
plot(t, Ptank*1e-3, "k", "LineWidth",1.5);
xlabel("Time [sec]");
ylabel("Actuator Stagnation Pressure [kPa]");
hold off

figure(2)
hold on
plot(t, M, "k", "LineWidth",1.5);
xlabel("Time [sec]");
ylabel("Mach Number");
hold off

figure(3)
hold on
plot(t, volume, "k", "LineWidth",1.5);
xlabel("Time [sec]");
ylabel("Volume [$m^3$]","Interpreter","latex");
hold off

figure(4)
hold on
plot(t, m, "k", "LineWidth",1.5);
xlabel("Time [sec]");
ylabel("mass [$kg$]","Interpreter","latex");
hold off

radius = r;
figure(6)
hold on
plot(radius, Ptank*1e-3, "k", "LineWidth",1.5);
xlabel("radius [m]");
ylabel("Actuator Stagnation Pressure [kPa]");
hold off

%Linear regression to fit a 2nd order polynomal to graph 6
samplePoints = radius(1):0.001:radius(end);
pInterp = interp1(radius,Ptank,samplePoints);
figure(7)
hold on
plot(pInterp,samplePoints)
theta = linearRegression(samplePoints.',pInterp.');

for i = 1:1:length(samplePoints)
    curveFitP(i) = theta(1)*(pInterp(i)^2) + theta(2)*pInterp(i) + theta(3);
end

plot(pInterp,curveFitP)
hold off


exp = [0	37.5e-3; 0.1	42.5e-3; 0.2	43.5e-3; 0.4	46e-3; 0.5	48e-3; 0.6	48.5e-3; 0.7	49e-3; 0.8	50.5e-3];

figure(8)
hold on
plot(t, radius, "k", "LineWidth",1.5);
plot(exp(:,1),exp(:,2), "k--", "LineWidth",1.5)
xlabel("Time [sec]","Interpreter","latex",fontsize=14);
ylabel("Radius [$m$]","Interpreter","latex",fontsize=14);
legend("Curve Fit", "Real Data","Interpreter","latex","Location", "southeast",fontsize=16);
hold off



%{
Model of the soft actuator's volume as a function of pressure
%}
function [v r] = updateVolume(P, l0, r0, t0, E)
    E1 = 58460;
    E2 = 10540;

    Pgage = P - 101325;

    a = E2*t0/(r0*(r0 + t0));
    b = E1*log((r0 + t0)/r0) - Pgage;
     c = -r0*Pgage;
 
    delta_r = (-b + sqrt(b.^2 - 4*a*c)) / (2 * a);
    r = r0 + delta_r;

    v = pi*(r)^2*l0;
end

%{
performs a second order curve fit by doing linear regression
%}
function theta = linearRegression(l,p)
    %with linear regressio
    phi = [p.^2 p ones(length(p),1)];
    
    theta = inv(phi.' * phi)*phi.' * l;
end
