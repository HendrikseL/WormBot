%{
Calculates the speed of delfation for the center soft actuator.
%}

%need to find how much pressure will change
clear
close all
tic

%Fluid Properties
Ptank = 101325;%starts at atmo pressure 
R = 287; %J/kgK
Pcomp = 10000 + Ptank; %working pressure + atmospheric
Tstag = 298; %K, isothermal model Tstag will not change
rhoTank = Ptank/(Tstag*R);
gamma = 1.4;

%geometric properties
r_init = (75/2)*1e-3; %m
l_init = 0.15; %initial length
thickness = 5e-3; %m
E = 69e3; %young's modulus, Pa
Avalve = 3.1416e-06; %area of the valve, circle of radius 10mm

%Assuming initially the solenoid is choked
M = real(sqrt( (2/(gamma-1)) * ( (Pcomp/Ptank)^((gamma-1)/gamma) -1 )));

if M > 1
    M =1;
end

%stagnation density of fluid in the compressor
rhoComp= ((Pcomp/Ptank)^(1/gamma))*rhoTank;
%inital volume of fully inflated tank
volume = updateVolume(Ptank, l_init,r_init,thickness, E); %m^3 

%initial mass
m = rhoTank * volume;

%under relaxation factor
a = 0.1;

%calculate Constant mass flow
T = Tstag * (1+ ((gamma-1)/2)*M^2)^(-1);
P = Pcomp * (1+ ((gamma-1)/2)*M^2)^((-gamma)/(gamma-1));
rho = P/(R*T);

%Mass flow
v = M * sqrt(gamma * R * T);
mdot = rho*v*Avalve;

deltat = 0.0001;
t = 0;
i = 1;
%main loop
while Ptank < Pcomp
   
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
        vol_update = updateVolume(Ptank(i+1), l_init,r_init,thickness, E);
        vol_inter(j+1) = (1-a)*vol_inter(j) + a*vol_update;

        rhoTank = m(i+1) / vol_inter(j+1);

        %convergence check
        if abs(vol_inter(j+1)-vol_inter(j)) < 1e-9
            flag = 0;
            volume(i+1) = vol_inter(j+1);
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

%{
Model of the soft actuator's volume as a function of pressure
%}
function v = updateVolume(P, l0, r0, t0, E)
    E1 = 58460;
    E2 = 10540;

    Pgage = P - 101325;
    c = r0^2/((r0 + t0)^2 - r0^2);

    sigma = (-E1 + sqrt(E1^2 + 4*E2*c*Pgage))/(2*E2);
    l = l0 + l0*sigma;
    v = pi*(r0)^2*l;

    % sigma = 1 + 2*r0^2/(E*(2*r0*t0 + t0^2))*Pgage;
    % l = l0*sigma;
    % v = pi*((r0)^2)*l;
end
