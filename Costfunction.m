function Cost = Costfunction(Np, Tim_step, X0, u, Vehicle_Type, Q, Xdes, R, Udes, F, Xa, G, Xnfa,Xnba)
%UNTITLED3 Summary of this function goes here
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
    
    Xp = [Pp,Vp];      % Predictive State
    
    Cost = 0;
    for i = 1:Np        %% 注意范数的定义问题， X'Q'QX
        Cost = Cost + (Xp(i,:)-Xdes(i,:))*Q*(Xp(i,:)-Xdes(i,:))' + ...
                (u(i)-Udes(i))*R*(u(i)-Udes(i)) + (Xp(i,:)-Xa(i,:))*F*(Xp(i,:)-Xa(i,:))'+ ...
                (Xp(i,:)-Xnfa(i,:))*G*(Xp(i,:)-Xnfa(i,:))'+(Xp(i,:)-Xnba(i,:))*G*(Xp(i,:)-Xnba(i,:))';               
    end
   
end

