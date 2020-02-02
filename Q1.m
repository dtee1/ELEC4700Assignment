% David Talson
% ELEC4700 Assignment 1

clear all
clc

global C

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.m_eff = 0.26*C.m_0;               % Effective mass of electrons 

regionWidth = 2e-7;                 % Nominal width of the region
regionLength = 1e-7;                % Nominal length of the region
T = 300;                            % Assumed temperature

vth = sqrt((2*C.kb*T)/C.m_eff);   % Thermal velocity as mean of magnitude of velocity
tmn = 2e-13;                        % Mean time between collisions
freepath = vth*2e-13;

numElectrons = 10;
electronXpos = rand(1,numElectrons).*regionWidth;
electronYpos = rand(1,numElectrons).*regionLength;
angle = rand(1,numElectrons).*2*pi;

% A negative value for the electron on the x-axis means the electron is
% moving to the left
% A negative value for the electron on the y-axis means the electron is
% moving downwards
electronXvel = vth.*cos(angle);
electronYvel = vth.*sin(angle);
deltaT = 1e-9/vth;
absTol = 0.0001e-7;
for i=1:1000
        for j = 1:numElectrons      
            if electronYpos(1,j) + electronYvel(1,j).*deltaT  >= regionLength || electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0
               electronYvel(1,j) = electronYvel(1,j)*-1;  
            end
            yposNew(i,j) = electronYpos(1,j) + electronYvel(1,j).*deltaT;
            if electronXpos(1,j) + electronXvel(1,j)*deltaT  >= regionWidth || electronXpos(1,j) + electronXvel(1,j)*deltaT <= 0
                if electronXpos(1,j) + electronXvel(1,j)*deltaT >= regionWidth
                    xposNew(i,j) = 0;
                else
                    xposNew(i,j) = regionWidth-0.01e-7;
                end
            else
                xposNew(i,j) = electronXpos(1,j) + electronXvel(1,j).*deltaT;
            end
        end
    
    electronXpos = xposNew(i,:);
    electronYpos = yposNew(i,:);
    
    % To calculate and plot the semiconductor temperature
    t1 = sqrt(electronXvel(1,:).^2 + electronYvel(1,:).^2);
    temperature(i) = ((mean(t1)^2)*C.m_eff)/(2*C.kb);

    
 
end

% Plot of the trajectories of the electrons
figure (1)
plot(xposNew,yposNew,'-','LineWidth',2)
xlim([-0.1e-7 2.1e-7])
ylim([-0.1e-7 1.1e-7])
grid on
xlabel('Electrons x position (m)')
ylabel('Electrons y position (m)')
title('A plot of path of 50 electrons in random motion')
figure(2)
plot(temperature)
grid on
xlabel('time')
ylabel('Temperature (K)')
title('A plot of temperature of semiconductor over time')

