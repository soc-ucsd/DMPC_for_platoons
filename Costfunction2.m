function Cost = Costfunction2(Np, Tim_step, X0, u, Vehicle_Type, Q, Xdes, R, F, Xa, G, Xnba)
% Cost function

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
    
    Udes = Radius/Eta*(Ca*Vp.^2 + Mass*g*f);
    U0 = Radius/Eta*(Ca*X0(2).^2 + Mass*g*f);
    
    Cost = (X0(1:2)-Xdes(1,:))*Q*(X0(1:2)-Xdes(1,:))' + ...
                (u(1)-U0)*R*(u(1)-U0) + (X0(1:2)-Xa(1,:))*F*(X0(1:2)-Xa(1,:))'+ ...
                (X0(1:2)-Xnba(1,:))*G*(X0(1:2)-Xnba(1,:))';                               % 第一步的优化值
    for i = 1:Np-1        %% 注意范数的定义问题， X'Q'QX
        Cost = Cost + (Xp(i,:)-Xdes(i+1,:))*Q*(Xp(i,:)-Xdes(i+1,:))' + ...
                (u(i+1)-Udes(i))*R*(u(i+1)-Udes(i)) + (Xp(i,:)-Xa(i+1,:))*F*(Xp(i,:)-Xa(i+1,:))'+ ...
                (Xp(i,:)-Xnba(i+1,:))*G*(Xp(i,:)-Xnba(i+1,:))';               
    end
   
end

