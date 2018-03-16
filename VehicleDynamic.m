function [PositionN, VelocityN, TorqueN] = VehicleDynamic(u,Tim_step,Position,Velocity,Torque,Mass,R,g,f,Eta,Ca,Tao)
% Vehicle dynamics

    PositionN = Position + Velocity*Tim_step;
    VelocityN = Velocity + 1/Mass *(Eta*Torque/R - Ca*Velocity^2 - Mass*g*f)*Tim_step;
    TorqueN = Torque - 1/Tao*Torque*Tim_step + 1/Tao*u*Tim_step;
    
end

