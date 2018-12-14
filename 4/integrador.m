function [xe1kk,P1kk] = integrador(xekk,Pkk,uk,Q,h)
    
    F = [-1 1; -0.2*xekk(1) 0];

    xp = [-xekk(1)+xekk(2); -0.1*xekk(1)^2 - 1 + uk];
    Pp = F*Pkk + Pkk*F' + Q;
    k1x = h*xp;
    k1p = h*Pp;
    
    xp = [-xekk(1)- k1x(1)/2 + xekk(2) + k1x(2)/2; -0.1*(xekk(1)+k1x(1)/2)^2 - 1 + uk];
    Pp = F*(Pkk + k1p/2) + (Pkk*k1p/2)*F' + Q;
    k2x = h*xp;
    k2p = h*Pp;
    
    xp = [-xekk(1)- k2x(1)/2 + xekk(2) + k2x(2)/2; -0.1*(xekk(1)+k2x(1)/2)^2 - 1 + uk];
    Pp = F*(Pkk + k2p/2) + (Pkk*k2p/2)*F' + Q;
    k3x = h*xp;
    k3p = h*Pp;
    
    xp = [-xekk(1)- k3x(1) + xekk(2) + k3x(2); -0.1*(xekk(1)+k3x(1))^2 - 1 + uk];
    Pp = F*(Pkk + k3p) + (Pkk*k3p)*F' + Q;
    k4x = h*xp;
    k4p = h*Pp;

    xe1kk = xekk + k1x/6 + k2x/3 + k2x/3 + k4x/6;
    P1kk = Pkk + k1p/6 + k2p/3 + k2p/3 + k4p/6;
end

