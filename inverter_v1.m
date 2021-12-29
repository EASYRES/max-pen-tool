function [ Sa, Sb, Sc, status] = inverter_v1( Va, Vb, Vc, Pref, Sa_old, Sb_old, Sc_old, Prated, ctrl_type)
% Va, Vb, Vc: phase voltages in complex form
% Pref: Active power reference (in kW)
% Sa_old, Sb_old, Sc_old: Inverter's complex powers of previous iteration
% (in kVA)

% Calculation of positive sequence voltage (in pu)
V1 = (Va+Vb*exp(1i*2*pi/3)+Vc*exp(1i*4*pi/3))/3;
v1 = abs(V1)/(400/sqrt(3));
k = 6; %integral term
Vt1 = 1.02; % Q(V) droop voltage threshold
Vt2 = 1.04; % Q(V) droop voltage threshold

switch ctrl_type
    case 1  % UPF
        Sa = Pref/3 + 1i*0;
        Sb = Pref/3 + 1i*0;
        Sc = Pref/3 + 1i*0;
        status = 0;
  
        
    case 2  % Q(P)
        if Pref < 0.5*Prated
            Sa = Pref/3 + 1i*0;
            Sb = Pref/3 + 1i*0;
            Sc = Pref/3 + 1i*0;
        elseif Pref >= 0.5*Prated && Pref < 0.9*Prated
            Sa = Pref/3 - 1i*(1.1*(Pref/3-0.5*Prated/3));
            Sb = Pref/3 - 1i*(1.1*(Pref/3-0.5*Prated/3));
            Sc = Pref/3 - 1i*(1.1*(Pref/3-0.5*Prated/3)); 
        else
            Sa = Pref/3 - 1i*0.44/3;
            Sb = Pref/3 - 1i*0.44/3;
            Sc = Pref/3 - 1i*0.44/3; 
        end
        status = 0;
            
    case 3
        slope_factor = -7.33;
        if (v1 <= Vt1)
            Qtemp = 0;
        elseif (Vt1 < v1 && v1 <= Vt2)
            Qtemp = slope_factor*(v1-Vt1)*Prated;
        else
            Qtemp = -0.44*Prated;
        end
   
        Sa_temp = Pref/3 + 1i*Qtemp/3;
        Sb_temp = Pref/3 + 1i*Qtemp/3;
        Sc_temp = Pref/3 + 1i*Qtemp/3;
        
        % Revised complex power output with the insertion of integral term k
        Sa = Sa_old + (Sa_temp-Sa_old)/k;
        Sb = Sb_old + (Sb_temp-Sb_old)/k;
        Sc = Sc_old + (Sc_temp-Sc_old)/k;
        status = 1;
    case 4
        
        

end


end

