%{
Calculates the speed of delfation for the center soft actuator.
%}

%need to find how much pressure will change
clear
close all
tic



%Fluid Properties
Patm = 99610; 
R = 287; %J/kgK
Pstag = 105859; %working pressure + atmospheric
Tstag = 298; %K, isothermal model Tstag will not change
rhoatm = Patm/(Tstag*R);
gamma = 1.4;

%geometric properties
r_init = (75/2)*1e-3; %m
l_init = 0.2; %initial length
thickness = 5e-3; %m
E = 69e3; %young's modulus, Pa
% Avalve = 3.1416e-06; %area of the valve, circle of radius 1mm
Avalve = 3.0660e-05;

%Assuming initially the solenoid is choked
M = real(sqrt( (2/(gamma-1)) * ( (Pstag/Patm)^((gamma-1)/gamma) -1 )));

if M > 1
    M =1;
end


%stagnation density of fluid in the tank
rhostag= ((Pstag/Patm)^(1/gamma))*rhoatm;
%inital volume of fully inflated tank
volume = updateVolume(Pstag, l_init,r_init,thickness, E, []); %m^3 

%initial mass
m = rhostag * volume;

%under relaxation factor
a = 0.1;

deltat = 0.0001;
t = 0;
i = 1;
%main loop
while Pstag > Patm
    %calculate static values
    T = Tstag * (1+ ((gamma-1)/2)*M(i)^2)^(-1);
    P = Pstag(i) * (1+ ((gamma-1)/2)*M(i)^2)^((-gamma)/(gamma-1));
    rho = P/(R*T);

    %Mass flow
    v(i) = M(i) * sqrt(gamma * R * T);
    mdot = rho*v(i)*Avalve;
   
    %update stagnation density and total mass
    m(i+1) = m(i) - (mdot * deltat);
    rhostag = m(i+1) / volume(i);
    
   
    %inner loop, checking for volume convergence
    flag = 0;
    j = 1;
    vol_inter = zeros(1000,1);
    vol_inter(1) = volume(i);
    while flag == 0 && j < 3000
        % update stagnation pressure
        Pstag(i+1) = Patm*(rhostag/rhoatm)^1.4;
    
        %update volume
        vol_update = updateVolume(Pstag(i+1), l_init,r_init,thickness, E, volume);
        vol_inter(j+1) = (1-a)*vol_inter(j) + a*vol_update;

        rhostag = m(i+1) / vol_inter(j+1);

        %convergence check
        if abs(vol_inter(j+1)-vol_inter(j)) < 1e-9
            flag = 1;
            volume(i+1) = vol_inter(j+1);
        end

        %Convergence iteration
        j = j +1;
    end

    %Update mach number
    M(i+1) = real(sqrt( (2/(gamma-1)) * ( (Pstag(i+1)/Patm)^((gamma-1)/gamma) -1 )));

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
plot(t, Pstag*1e-3, "k", "LineWidth",2);
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

figure(5)
hold on
plot(volume, Pstag*1e-3, "k", "LineWidth",1.5);
xlabel("Volume [m^3]");
ylabel("Actuator Stagnation Pressure [kPa]");
hold off

lengths = (volume/(pi*(r_init)^2));
figure(6)
hold on
plot(lengths, Pstag*1e-3, "k", "LineWidth",1.5);
xlabel("Length [m]");
ylabel("Actuator Stagnation Pressure [kPa]");
hold off

%Linear regression to fit a 2nd order polynomal to graph 6
samplePoints = lengths(end):0.001:lengths(1);
pInterp = interp1(lengths,Pstag,samplePoints);
figure(7)
hold on
plot(pInterp,samplePoints)
theta = linearRegression(samplePoints.',pInterp.');

for i = 1:1:length(samplePoints)
    curveFitP(i) = theta(1)*(pInterp(i)^2) + theta(2)*pInterp(i) + theta(3);
end

plot(pInterp,curveFitP)
hold off



modelledPressure(:,1) = t.';
modelledPressure(:,2) = Pstag.';

%read in data
file = "yank_5.csv";
fid = fopen(file);
tline = fgetl(fid);

%Read data
i=1;
while tline ~= -1

    if ~strcmp(tline,"b'Valves set.'")
        c = textscan(tline,'%f %f %f %f %f %f %f %f','Delimiter',',');
        c = string(c);
    
        %points
        try
            pressure(i,2) = round(sscanf(sprintf(' %s',c{1,2}),'%f',[1,Inf]),9);
        
            if pressure(i,2) < 115000
                pressure(i,1) = round(sscanf(sprintf(' %s',c{1,1}),'%f',[1,Inf]),9)/(1e6);
                i = i+1;
            end
        catch
            %do nothing
        end
    end

    tline = fgetl(fid);
    
end
fclose(fid);

% load("yank5_pressure_dynamic.mat");
modelledPressure(end+1,:) = [3 modelledPressure(length(modelledPressure)-1,2)];

figure(8)
hold on
plot(pressure(:,1)-10.92,pressure(:,2)/1000,'k',LineWidth=1.25)
plot(modelledPressure(:,1),modelledPressure(:,2)/1000,'k--',LineWidth=1.25)
xlim([0 1])
xlabel("Time [sec]",fontsize=14)
ylabel("Stagnation Pressure [kPa]",fontsize=14)
legend("Real Data Curve", "Simulated Curve",fontsize=16)
hold off

%{
Model of the soft actuator's volume as a function of pressure
%}
function v = updateVolume(P, l0, r0, t0, E, volume)
    E1 = 58460;
    E2 = 10540;
    
    Pgage = P - 101325;
    c = r0^2/((r0 + t0)^2 - r0^2);

    sigma = (-E1 + sqrt(E1^2 + 4*E2*c*Pgage))/(2*E2);
    l = l0 + l0*sigma; %elong;



    % elong = 3.6271e-10*Pgage^2  + 2.1408e-06*Pgage  + 3.7846e-04;
    % l = l0 + elong;
    v = pi*(r0)^2*l;

    % sigma = 1 + 2*r0^2/(E*(2*r0*t0 + t0^2))*Pgage;
    % l = l0*sigma;
    % v = pi*((r0)^2)*l;
end


%{
performs a second order curve fit by doing linear regression
%}
function theta = linearRegression(l,p)
    %with linear regressio
    phi = [p.^2 p ones(length(p),1)];
    
    theta = inv(phi.' * phi)*phi.' * l;
end


%% System dynamics with no friction (state space form)
% Input should be a vector with x1, v1, x2, v2
% l0, r0, t0 are geometry of the unrestrained portion of the tube
% stress_func returns a stress given a strain
function xdot = dynamics_ss_fixed(x, P, l0, r0, t0, m1, m2, mc, b, stress_func)
  V0 = l0*pi*((r0 + t0)^2 - r0^2);
  strain = (x(2) - l0) / l0;
  xdot = [x(2);
          (-b*x(2) - V0/l0*stress_func(strain))/(m2 + 1/6*mc)] + ...
          pi*r0^2*P*[0 1/(m2 + 1/6*mc)]';
  return;
end

function stress = stress_func(strain)
    stress = 62150*strain^0.7036;
end