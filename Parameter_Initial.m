
clc;clear;close all;

%% Parameter Initialation

Num_veh  = 7;                            % The number of vehicles in a platoon; 
Tim_step = 0.1;                          % Time step;
Time_sim = 10;                           % Time length for simulation
Num_step = Time_sim/Tim_step;            % Simulation setps

Mass = 1000 + zeros(Num_veh,1) + 1000*rand(Num_veh,1);  % Vehilce mass
Tao  = 0.5 +(Mass - 1000)/1000 * 0.3;   % Time lag
f    = 0.01;                            % rolling friction
Eta  = 0.96;                             % Efficency
g    = 9.8;                             % gravitity 
Ca   = 0.98 + (Mass - 1000)/1000 * 0.2;                  % ·ç×è 
% Ca   = 1/2*0.4*2*1.23;                  % ·ç×è 
R    = 0.3 +(Mass - 2000)/2000 * 0.1;   % °ë¾¶

%% Acceleration bounds -6,6
AccMax = 6; AccMin = -6;
Torquebound = zeros(6,2);               %   [low up]
for i = 1:Num_veh
    Torquebound(i,1) = Mass(i)*AccMin*R(i)/Eta;
    Torquebound(i,2) = Mass(i)*AccMax*R(i)/Eta;
end
