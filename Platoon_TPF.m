
%% Code for the paper
% Title:   Distributed model predictive control for heterogeneous vehicle platoons under unidirectional topologies
% Authors: Zheng, Yang, Shengbo Eben Li, Keqiang Li, Francesco Borrelli, and J. Karl Hedrick. 
% Journal: IEEE Transactions on Control Systems Technology 25, no. 3 (2017): 899-910.


%% DMPC for platoons with TPF topology
clc;clear;close all;
load PlatoonParameter.mat  % This set of parameters were used in the paper


%% Initial Virables 
Postion  = zeros(Num_step,Num_veh);     % postion of each vehicle;
Velocity = zeros(Num_step,Num_veh);     % velocity of each vehicle;
Torque   = zeros(Num_step,Num_veh);     % Braking or Tracking Torque of each vehicle;
U        = zeros(Num_step,Num_veh);     % Desired Braking or Tracking Torque of each vehicle;

Cost    = zeros(Num_step,Num_veh);      % Cost function
Exitflg = zeros(Num_step,Num_veh);      % Stop flag - solvers

% Leading vehicle
d  = 20;                                % Desired spacing
a0 = zeros(Num_step,1); 
v0 = zeros(Num_step,1); 
x0 = zeros(Num_step,1);


% Transient process of leader, which is given in advance
v0(1) = 20; a0(1/Tim_step+1:2/Tim_step) = 2; 
for i = 2:Num_step
    v0(i) = v0(i-1)+a0(i)*Tim_step; 
    x0(i) = x0(i-1)+v0(i)*Tim_step;    
end
 
% Zero initial error for the followers
for i = 1:Num_veh
    Postion(1,i) = x0(1)-i*d;
    Velocity(1,i) = 20;             
    Torque(1,i) = (Mass(i)*g*f + Ca(i)*Velocity(1,i)^2)*R(i)/Eta;
end



%% Iterative Simulation 
% TPF topology --> Fi > Gi+1
% Q1 : leader weighted matrix for state; 
% R1 --> leader weighted matrix for control input
% Fi --> 自身偏离的惩罚
% Gi --> 偏离邻域车辆状态的惩罚

%  MPC weighted matrix initial
F1 = 10*eye(2); G1 = 0;          Q1 = 10*eye(2);R1 = 1;     
F2 = 10*eye(2); G2 = 10/2*eye(2);Q2 = 10*eye(2); R2 = 1;     
F3 = 10*eye(2); G3 = 10/2*eye(2);Q3 = 0*eye(2); R3 = 1;      
F4 = 10*eye(2); G4 = 10/2*eye(2);Q4 = 0*eye(2); R4 = 1;     
F5 = 10*eye(2); G5 = 10/2*eye(2);Q5 = 0*eye(2); R5 = 1;     
F6 = 10*eye(2); G6 = 10/2*eye(2);Q6 = 0*eye(2); R6 = 1;      
F7 = 10*eye(2); G7 = 10/2*eye(2);Q7 = 0*eye(2); R7 = 1;       

% Distributed MPC assumed state
Np = 20;                      % 预测步长
Pa = zeros(Np,Num_veh);       % Assumed postion of each vehicle;
Va = zeros(Np,Num_veh);       % Assumed velocity of each vehicle;
ua = zeros(Np,Num_veh);       % Assumed  Braking or Tracking Torque input of each vehicle;

Pa_next = zeros(Np+1,Num_veh);  % 1（0）：为上一时刻的状态Assumed postion of each vehicle at the newt time step;
Va_next = zeros(Np+1,Num_veh);  % Assumed velocity of each vehicle at the newt time step;
ua_next = zeros(Np+1,Num_veh);  % Assumed Braking or Tracking Torque of each vehicle at the newt time step;

% Initialzie the assumed state for the first computation: constant speed
for i = 1:Num_veh
    ua(:,i) = Torque(1,i);
    Pa(1,i) = Postion(1,i);                % 假设的第一个点  为文章中的 k=0处，为当前车辆的状态；
    Va(1,i) = Velocity(1,i);
    Ta(1,i) = Torque(1,i);
    for j = 1:Np
        [Pa(j+1,i),Va(j+1,i),Ta(j+1,i)] = VehicleDynamic(ua(j,i),Tim_step,Pa(j,i),Va(j,i),Ta(j,i),Mass(i),R(i),g,f,Eta,Ca(i),Tao(i));
    end    
end

tol_opt = 1e-5;
options = optimset('Display','off','TolFun', tol_opt, 'MaxIter', 2000,...
                'LargeScale', 'off', 'RelLineSrchBnd', [], 'RelLineSrchBndDuration', 1);
 
 %%  调试用
 % 终端状态 
 Xend = zeros(Num_step,Num_veh); Vend = zeros(Num_step,Num_veh);
 
 %%  循环仿真

for i = 2:Num_step - Np
    
    fprintf('\n Steps i= %d\n',i)
    
    % Solve optimization problem
    tic
    %% Vehicle one
    Vehicle_Type = [Mass(1),R(1),g,f,Eta,Ca(1),Tao(1)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,1),Velocity(i-1,1),Torque(i-1,1)];                % the vehicle variable in the last time
    Pd = x0(i-1:i+Np-1) - d;  Vd = v0(i-1:i+Np-1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,1),Va(:,1)];                                             % 自己预期的行为，传递给下一辆车
    Xnba = zeros(Np+1,2);                                               % 1：为上一时刻的状态
   
    u0 = ua(:,1);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(1,1)*ones(Np,1); ub = Torquebound(1,2)*ones(Np,1);      % 控制量上下界              
    Pnp = Pd(end,1); Vnp = Vd(end,1);   % 终端约束
    Xend(i,1) = Pnp; Vend(i,1) = Vnp; Tnp = (Ca(1)*Vnp.^2 + Mass(1)*g*f)/Eta*R(1);
    % MPC 优化求解
    [u, Cost(i,1), Exitflg(i,1), output] = fmincon(@(u) Costfunction2( Np, Tim_step, X0 ,u, Vehicle_Type,Q1,Xdes,R1,F1,Xa,G1,Xnba), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,1) = u(1);
    [Postion(i,1),Velocity(i,1),Torque(i,1)] = VehicleDynamic(U(i,1),Tim_step,Postion(i-1,1),Velocity(i-1,1),Torque(i-1,1),Mass(1),R(1),g,f,Eta,Ca(1),Tao(1));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,1),Velocity(i,1),Torque(i,1)];   
    ua(1:Np-1,1) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,1),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(1),R(1),g,f,Eta,Ca(1),Tao(1));
    end    
    ua(Np,1) = (Ca(1)*Temp(Np,2).^2 + Mass(1)*g*f)/Eta*R(1);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,1),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(1),R(1),g,f,Eta,Ca(1),Tao(1));
    Pa_next(:,1) = Temp(:,1);
    Va_next(:,1) = Temp(:,2);
    toc
    
    
    %% Vehicle two
    tic
    Vehicle_Type = [Mass(2),R(2),g,f,Eta,Ca(2),Tao(2)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,2),Velocity(i-1,2),Torque(i-1,2)];                % the vehicle variable in the last time
    Pd = x0(i-1:i+Np-1) - 2*d;  Vd = v0(i-1:i+Np-1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,2),Va(:,2)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,1) - d, Va(:,1)];                                               % 1：为上一时刻的状态
   
    u0 = ua(:,2);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(2,1)*ones(Np,1); ub = Torquebound(2,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Pd(end))/2; Vnp = (Xnfa(end,2)+Vd(end))/2;   % 终端约束
    Xend(i,2) = Pnp; Vend(i,2) = Vnp; Tnp = (Ca(2)*Vnp.^2 + Mass(2)*g*f)/Eta*R(2);
    % MPC 优化求解
    [u, Cost(i,2), Exitflg(i,2), output] = fmincon(@(u) Costfunction2( Np, Tim_step, X0 ,u, Vehicle_Type,Q2,Xdes,R2,F2,Xa,G2,Xnfa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,2) = u(1);
    [Postion(i,2),Velocity(i,2),Torque(i,2)] = VehicleDynamic(U(i,2),Tim_step,Postion(i-1,2),Velocity(i-1,2),Torque(i-1,2),Mass(2),R(2),g,f,Eta,Ca(2),Tao(2));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,2),Velocity(i,2),Torque(i,2)]; 
    ua(1:Np-1,2) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,2),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(2),R(2),g,f,Eta,Ca(2),Tao(2));
    end    
    
    ua(Np,2) = (Ca(2)*Temp(Np,2).^2 + Mass(2)*g*f)/Eta*R(2);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,2),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(2),R(2),g,f,Eta,Ca(2),Tao(2));
    Pa_next(:,2) = Temp(:,1);
    Va_next(:,2) = Temp(:,2);
    toc
    
    
    
    %% vehicle three
    tic
    Vehicle_Type = [Mass(3),R(3),g,f,Eta,Ca(3),Tao(3)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,3),Velocity(i-1,3),Torque(i-1,3)];                % the vehicle variable in the last time
    Pd = zeros(Np+1,1);  Vd = zeros(Np+1,1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,3),Va(:,3)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,2) - d, Va(:,2)];                                               % 1：为上一时刻的状态
    Xnffa = [Pa(:,1) - 2*d, Va(:,1)]; 
    
    u0 = ua(:,3);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(3,1)*ones(Np,1); ub = Torquebound(3,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Xnffa(end,1))/2; Vnp = (Xnfa(end,2)+Xnffa(end,2))/2;   % 终端约束
    Xend(i,3) = Pnp; Vend(i,3) = Vnp; Tnp = (Ca(3)*Vnp.^2 + Mass(3)*g*f)/Eta*R(3);
    % MPC 优化求解
    [u, Cost(i,3), Exitflg(i,3), output] = fmincon(@(u) Costfunction1( Np, Tim_step, X0 ,u, Vehicle_Type,Q3,Xdes,R3,F3,Xa,G3,Xnfa,Xnffa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,3) = u(1);
    [Postion(i,3),Velocity(i,3),Torque(i,3)] = VehicleDynamic(U(i,3),Tim_step,Postion(i-1,3),Velocity(i-1,3),Torque(i-1,3),Mass(3),R(3),g,f,Eta,Ca(3),Tao(3));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,3),Velocity(i,3),Torque(i,3)];    
    ua(1:Np-1,3) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,3),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(3),R(3),g,f,Eta,Ca(3),Tao(3));
    end   
    
    ua(Np,3) = (Ca(3)*Temp(Np,2).^2 + Mass(3)*g*f)/Eta*R(3);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,3),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(3),R(3),g,f,Eta,Ca(3),Tao(3));
    Pa_next(:,3) = Temp(:,1);
    Va_next(:,3) = Temp(:,2);
    toc
    
    
    
    %% vehicle four
    tic
    Vehicle_Type = [Mass(4),R(4),g,f,Eta,Ca(4),Tao(4)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,4),Velocity(i-1,4),Torque(i-1,4)];                % the vehicle variable in the last time
    Pd = zeros(Np+1,1);  Vd = zeros(Np+1,1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,4),Va(:,4)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,3) - d, Va(:,3)];                                               % 1：为上一时刻的状态
    Xnffa = [Pa(:,2) - 2*d, Va(:,2)];
   
    u0 = ua(:,4);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(4,1)*ones(Np,1); ub = Torquebound(4,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Xnffa(end,1))/2; Vnp = (Xnfa(end,2)+Xnffa(end,2))/2;   % 终端约束
    Xend(i,4) = Pnp; Vend(i,4) = Vnp; Tnp = (Ca(4)*Vnp.^2 + Mass(4)*g*f)/Eta*R(4);
    % MPC 优化求解
    [u, Cost(i,4), Exitflg(i,4), output] = fmincon(@(u) Costfunction1( Np, Tim_step, X0 ,u, Vehicle_Type,Q3,Xdes,R3,F3,Xa,G3,Xnfa,Xnffa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,4) = u(1);
    [Postion(i,4),Velocity(i,4),Torque(i,4)] = VehicleDynamic(U(i,4),Tim_step,Postion(i-1,4),Velocity(i-1,4),Torque(i-1,4),Mass(4),R(4),g,f,Eta,Ca(4),Tao(4));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,4),Velocity(i,4),Torque(i,4)];  
    ua(1:Np-1,4) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,4),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(4),R(4),g,f,Eta,Ca(4),Tao(4));
    end
    
    ua(Np,4) = (Ca(4)*Temp(Np,2).^2 + Mass(4)*g*f)/Eta*R(4);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,4),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(4),R(4),g,f,Eta,Ca(4),Tao(4));
    Pa_next(:,4) = Temp(:,1);
    Va_next(:,4) = Temp(:,2);
    toc
    
    
    
      %% vehicle five
    tic
    Vehicle_Type = [Mass(5),R(5),g,f,Eta,Ca(5),Tao(5)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,5),Velocity(i-1,5),Torque(i-1,5)];                % the vehicle variable in the last time
    Pd = zeros(Np+1,1);  Vd = zeros(Np+1,1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,5),Va(:,5)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,4) - d, Va(:,4)];                                               % 1：为上一时刻的状态
    Xnffa = [Pa(:,3) - 2*d, Va(:,3)];
   
    u0 = ua(:,5);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(5,1)*ones(Np,1); ub = Torquebound(5,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Xnffa(end,1))/2; Vnp = (Xnfa(end,2)+Xnffa(end,2))/2;   % 终端约束
    Xend(i,5) = Pnp; Vend(i,5) = Vnp; Tnp = (Ca(5)*Vnp.^2 + Mass(5)*g*f)/Eta*R(5);
    % MPC 优化求解
    [u, Cost(i,5), Exitflg(i,5), output] = fmincon(@(u) Costfunction1( Np, Tim_step, X0 ,u, Vehicle_Type,Q3,Xdes,R3,F3,Xa,G3,Xnfa,Xnffa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,5) = u(1);
    [Postion(i,5),Velocity(i,5),Torque(i,5)] = VehicleDynamic(U(i,5),Tim_step,Postion(i-1,5),Velocity(i-1,5),Torque(i-1,5),Mass(5),R(5),g,f,Eta,Ca(5),Tao(5));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,5),Velocity(i,5),Torque(i,5)];   
    ua(1:Np-1,5) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,5),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(5),R(5),g,f,Eta,Ca(5),Tao(5));
    end
    
    ua(Np,5) = (Ca(5)*Temp(Np,2).^2 + Mass(5)*g*f)/Eta*R(5);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,5),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(5),R(5),g,f,Eta,Ca(5),Tao(5));
    Pa_next(:,5) = Temp(:,1);
    Va_next(:,5) = Temp(:,2);

    toc
    
    %% vehicle six
    tic
    Vehicle_Type = [Mass(6),R(6),g,f,Eta,Ca(6),Tao(6)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,6),Velocity(i-1,6),Torque(i-1,6)];                % the vehicle variable in the last time
    Pd = zeros(Np+1,1);  Vd = zeros(Np+1,1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,6),Va(:,6)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,5) - d, Va(:,5)];                                               % 1：为上一时刻的状态
    Xnffa = [Pa(:,4) - 2*d, Va(:,4)];
    
    u0 = ua(:,6);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(6,1)*ones(Np,1); ub = Torquebound(6,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Xnffa(end,1))/2; Vnp = (Xnfa(end,2)+Xnffa(end,2))/2;   % 终端约束
    Xend(i,6) = Pnp; Vend(i,6) = Vnp; Tnp = (Ca(6)*Vnp.^2 + Mass(6)*g*f)/Eta*R(6);
    % MPC 优化求解
    [u, Cost(i,6), Exitflg(i,6), output] = fmincon(@(u) Costfunction1( Np, Tim_step, X0 ,u, Vehicle_Type,Q3,Xdes,R3,F3,Xa,G3,Xnfa,Xnffa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,6) = u(1);
    [Postion(i,6),Velocity(i,6),Torque(i,6)] = VehicleDynamic(U(i,6),Tim_step,Postion(i-1,6),Velocity(i-1,6),Torque(i-1,6),Mass(6),R(6),g,f,Eta,Ca(6),Tao(6));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,6),Velocity(i,6),Torque(i,6)]; 
    ua(1:Np-1,6) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,6),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(6),R(6),g,f,Eta,Ca(6),Tao(6));
    end
     
    ua(Np,6) = (Ca(6)*Temp(Np,2).^2 + Mass(6)*g*f)/Eta*R(6);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,6),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(6),R(6),g,f,Eta,Ca(6),Tao(6));
    Pa_next(:,6) = Temp(:,1);
    Va_next(:,6) = Temp(:,2);

    toc
    
     %% vehicle seven
    tic
    Vehicle_Type = [Mass(7),R(7),g,f,Eta,Ca(7),Tao(7)];                 % the vehicle parameters ： Mass,R,g,f,Eta,Ca(i),Tao, 
    X0 = [Postion(i-1,7),Velocity(i-1,7),Torque(i-1,7)];                % the vehicle variable in the last time
    Pd = zeros(Np+1,1);  Vd = zeros(Np+1,1);                      % 共Np+1个点，注意下角标，i-1 代表上一时刻的状态， i代表当前需要优化求解的状态
    Xdes = [Pd,Vd];  % Udes = Td;                                       % 第一辆车的期望行为
    Xa = [Pa(:,7),Va(:,7)];                                             % 自己预期的行为，传递给下一辆车
    Xnfa = [Pa(:,6) - d, Va(:,6)];                                               % 1：为上一时刻的状态
    Xnffa = [Pa(:,5) - 2*d, Va(:,5)];
    
    u0 = ua(:,7);   % 起始搜索点    
    A = [];b = []; Aeq = []; beq = [];                                       % 没有线性约束
    lb = Torquebound(7,1)*ones(Np,1); ub = Torquebound(7,2)*ones(Np,1);      % 控制量上下界              
    Pnp = (Xnfa(end,1)+Xnffa(end,1))/2; Vnp = (Xnfa(end,2)+Xnffa(end,2))/2;   % 终端约束
    Xend(i,7) = Pnp; Vend(i,7) = Vnp; Tnp = (Ca(7)*Vnp.^2 + Mass(7)*g*f)/Eta*R(7);
    % MPC 优化求解
    [u, Cost(i,7), Exitflg(i,7), output] = fmincon(@(u) Costfunction1( Np, Tim_step, X0 ,u, Vehicle_Type,Q3,Xdes,R3,F3,Xa,G3,Xnfa,Xnffa), ...
        u0, A, b, Aeq, beq, lb, ub, @(u) Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,Pnp,Vnp,Tnp),options); 
    
    % 车辆往前走一步
    U(i,7) = u(1);
    [Postion(i,7),Velocity(i,7),Torque(i,7)] = VehicleDynamic(U(i,7),Tim_step,Postion(i-1,7),Velocity(i-1,7),Torque(i-1,7),Mass(7),R(7),g,f,Eta,Ca(7),Tao(7));
    
    % 这个地方需要注意,下一阶段的assumed state, 在t+1时刻预测Np自身的状态
    Temp = zeros(Np+1,3);
    Temp(1,:) = [Postion(i,7),Velocity(i,7),Torque(i,7)];
    ua(1:Np-1,7) = u(2:Np);
    for j = 1:Np-1
        [Temp(j+1,1),Temp(j+1,2),Temp(j+1,3)] = VehicleDynamic(ua(j,7),Tim_step,Temp(j,1),Temp(j,2),Temp(j,3),Mass(7),R(7),g,f,Eta,Ca(7),Tao(7));
    end

    ua(Np,7) = (Ca(7)*Temp(Np,2).^2 + Mass(7)*g*f)/Eta*R(7);
    [Temp(Np+1,1),Temp(Np+1,2),Temp(Np+1,3)] = VehicleDynamic(ua(Np,7),Tim_step,Temp(Np,1),Temp(Np,2),Temp(Np,3),Mass(7),R(7),g,f,Eta,Ca(7),Tao(7));
    Pa_next(:,7) = Temp(:,1);
    Va_next(:,7) = Temp(:,2);

    toc

    
    %% 跟新交换数据矩阵
    Pa = Pa_next;
    Va = Va_next;
       
end

FigurePlot