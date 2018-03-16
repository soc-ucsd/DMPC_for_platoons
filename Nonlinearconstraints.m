function [C, Ceq] = Nonlinearconstraints(Np, Tim_step, X0, u, Vehicle_Type,x0,v0,T0)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    Pp = zeros(Np,1);     % Predictive Position
    Vp = zeros(Np,1);     % Predictive Velocity
    Tp = zeros(Np,1);     % Predictive Torque
    
    Mass = Vehicle_Type(1);Radius = Vehicle_Type(2); g = Vehicle_Type(3);f = Vehicle_Type(4);
    Eta = Vehicle_Type(5);Ca = Vehicle_Type(6);Tao = Vehicle_Type(7);
    
    [Pp(1),Vp(1),Tp(1)] = VehicleDynamic(u(1),Tim_step,X0(1),X0(2),X0(3),Mass,Radius,g,f,Eta,Ca,Tao);
    for i = 1:Np-1
        [Pp(i+1),Vp(i+1),Tp(i+1)] = VehicleDynamic(u(i+1),Tim_step,Pp(i),Vp(i),Tp(i),Mass,Radius,g,f,Eta,Ca,Tao);
    end
    
    Xp = [Pp,Vp,Tp];      % Predictive State
    
    C = [];
    Ceq = [Pp(Np) - x0;Vp(Np)-v0; Tp(Np) - T0];%Radius*(Ca*v0^2 + Mass*g*f)/Eta];   


end

