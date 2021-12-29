function max_pen_calc( hashmap, num_pv, pwr_pv_init, dstr_fctr, idx, ctrl_type, ti )

volt = 1;
thermal = 1;
cnvrgc_fail = 1;


tap_set = 1; % Zero tap position

load([pwd,filesep,'Data',filesep,'Consumption and generation data'])
load([pwd,filesep,'Data',filesep,'Network Data'],'Busname','Cablename','Linename','PQLoadname','Sts','Cable','Line','Transf');

for i = 1 : num_pv
    PV_name(i) = {['PV',num2str(i)]};
end

Pmppt=dstr_fctr*pwr_pv_init*Pprofile*ones( 1 , num_pv );

% First power flow
% Initialization of the OpenDSS
DSSObj = actxserver('OpenDSSEngine.DSS');
Start = DSSObj.Start(0);
DSSCircuit = DSSObj.ActiveCircuit;
DSSText = DSSObj.Text;
DSSObj.ClearAll;
% path of the .dss file
DSSText.Command = ['compile ',pwd,filesep,'Data',filesep,'Examined_Network_',num2str(idx),'.dss'];

% Tap Determination
DSSText.Command = ['Transformer.Transf2ph1.Taps=[',num2str(tap_set),',1]'];
DSSText.Command = ['Transformer.Transf2ph2.Taps=[',num2str(tap_set),',1]'];
DSSText.Command = ['Transformer.Transf2ph3.Taps=[',num2str(tap_set),',1]'];


for i = 1 : length(PQLoadname)
    DSSCircuit.SetActiveElement(['load.',PQLoadname{i}]);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(PQLoad(ti,i)));
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(PQLoad(ti,i)));
end

for i = 1 : length(PV_name)
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'a']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);
    
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'b']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);
    
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'c']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);

end

% First Load Flow Solution
DSSSolution=DSSCircuit.Solution;
DSSSolution.Solve;

% Voltage extraction
if DSSSolution.Converged==1
    Voltage_temp=voltage_extraction(Busname,DSSCircuit.AllNodename,DSSCircuit.AllBusVolts);
else
    display('Cannot converge');
     cnvrgc_fail = 0;
end


while volt && thermal && cnvrgc_fail
    
    pwr_pv_init = pwr_pv_init+1;
    Pmppt = dstr_fctr*pwr_pv_init*Pprofile*ones( 1 , num_pv );
    display(['Allocation number= ',num2str(idx),', Installed power = ',num2str(pwr_pv_init)]);
    
    for i = 1 : length(PQLoadname)
        DSSCircuit.SetActiveElement(['load.',PQLoadname{i}]);  
        DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(PQLoad(ti,i)));
        DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(PQLoad(ti,i)));
    end

    for i = 1 : length(PV_name)
        DSSCircuit.SetActiveElement(['generator.',PV_name{i},'a']);  
        DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
        DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);

        DSSCircuit.SetActiveElement(['generator.',PV_name{i},'b']);  
        DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
        DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);

        DSSCircuit.SetActiveElement(['generator.',PV_name{i},'c']);  
        DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(ti,i))/3);
        DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);
    end

    % First Load Flow Solution
    DSSSolution=DSSCircuit.Solution;
    DSSSolution.Solve;

    % Voltage extraction
    if DSSSolution.Converged==1
        Voltage_temp=voltage_extraction(Busname,DSSCircuit.AllNodename,DSSCircuit.AllBusVolts);
    else
        display('Cannot converge');
        cnvrgc_fail = 0;
    end

        
    %% Implementation of voltage control strategies
    cnvrgc = 1;
    iter = 1;
    status = 1;
    Pgen = Pmppt( ti, : );
    Qgen = zeros( 1, num_pv );
    S.a = Pgen/3+i*Qgen/3;
    S.b = Pgen/3+i*Qgen/3;
    S.c = Pgen/3+i*Qgen/3;
    
    
    while cnvrgc && iter<=Sts.localiter && status && cnvrgc_fail
        Stemp_a = zeros( 1, num_pv );
        Stemp_b = zeros( 1, num_pv );
        Stemp_c = zeros( 1, num_pv );

        for i = 1 : num_pv
            [ Sa, Sb, Sc, status ] = inverter_v1( Voltage_temp.VAN_vec(end,hashmap(i)), Voltage_temp.VBN_vec(end,hashmap(i)), Voltage_temp.VCN_vec(end,hashmap(i)), Pgen(i), S.a(end,i), S.b(end,i), S.c(end,i), dstr_fctr*pwr_pv_init, ctrl_type );
            Stemp_a(i) = Sa;
            Stemp_b(i) = Sb;
            Stemp_c(i) = Sc;

            DSSCircuit.SetActiveElement(['generator.',PV_name{i},'a']);  
            DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Sa));
            DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(Sa));

            DSSCircuit.SetActiveElement(['generator.',PV_name{i},'b']);
            DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Sb));
            DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(Sb));

            DSSCircuit.SetActiveElement(['generator.',PV_name{i},'c']);
            DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Sc));
            DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(Sc));
         end
        % Array that contains the injections of each inverter for all iterations        
        S.a=[S.a;Stemp_a];
        S.b=[S.b;Stemp_b];
        S.c=[S.c;Stemp_c];

         % Load flow Solution
        DSSSolution = DSSCircuit.Solution;
        DSSSolution.Solve;

         % Voltage extraction
        if DSSSolution.Converged == 1

            V_temp=voltage_extraction(Busname,DSSCircuit.AllNodename,DSSCircuit.AllBusVolts);

            Voltage_temp.VAN_vec=[Voltage_temp.VAN_vec;V_temp.VAN_vec];
            Voltage_temp.VBN_vec=[Voltage_temp.VBN_vec;V_temp.VBN_vec];
            Voltage_temp.VCN_vec=[Voltage_temp.VCN_vec;V_temp.VCN_vec];
            Voltage_temp.V1=V_temp.V1;
            Voltage_temp.V1_mag_pu=V_temp.V1_mag_pu;
            Voltage_temp.Vmax=V_temp.Vmax;
            Voltage_temp.Vmax_pu=V_temp.Vmax_pu;

        else

            display('OpenDSS cannot converge.');
            cnvrgc_fail = 0;
        end 
        SeqVoltage_temp.V0_vec=(Voltage_temp.VAN_vec+Voltage_temp.VBN_vec+Voltage_temp.VCN_vec)/3;
        SeqVoltage_temp.V1_vec=(Voltage_temp.VAN_vec+Voltage_temp.VBN_vec*exp(1i*2*pi/3)+Voltage_temp.VCN_vec*exp(1i*4*pi/3))/3;
        SeqVoltage_temp.V2_vec=(Voltage_temp.VAN_vec+Voltage_temp.VBN_vec*exp(1i*4*pi/3)+Voltage_temp.VCN_vec*exp(1i*2*pi/3))/3;

        % Check of the convergence criterion.
        summ=0;
        for bus=1:length(Busname)
            if abs(SeqVoltage_temp.V0_vec(iter+1,bus)-SeqVoltage_temp.V0_vec(iter,bus))<Sts.localtol && abs(SeqVoltage_temp.V1_vec(iter+1,bus)-SeqVoltage_temp.V1_vec(iter,bus))<Sts.localtol && abs(SeqVoltage_temp.V2_vec(iter+1,bus)-SeqVoltage_temp.V2_vec(iter,bus))<Sts.localtol
                summ=summ+1;
            else  
                break;
            end       
        end

        if summ==length(Busname)
            cnvrgc = 0;
        end

        iter=iter+1;

        if iter==Sts.localiter
            display(['Maximum iterations (',num2str(Sts.localiter),') reached. Local control cannot converge.']);
        end
    end
    
    
    %% Check for voltage or thermal violation
    cur_viol=0;

    for i=1:length(Linename)
        DSSCircuit.SetActiveElement(['line.',Linename{i}]);  
        curline1=DSSCircuit.ActiveElement.Currents;
        currentline1pha=curline1(1)+1i*curline1(2);
        currentline1phb=curline1(3)+1i*curline1(4);
        currentline1phc=curline1(5)+1i*curline1(6);
        currentline1phn=curline1(7)+1i*curline1(8);
        current.a_mang(i,1)=abs(currentline1pha);
        current.a_ang(i,1)=rad2deg(angle(currentline1pha));
        current.b_mang(i,1)=abs(currentline1phb);
        current.b_ang(i,1)=rad2deg(angle(currentline1phb));
        current.c_mang(i,1)=abs(currentline1phc);
        current.c_ang(i,1)=rad2deg(angle(currentline1phc));
        current.n_mang(i,1)=abs(currentline1phn);
        current.n_ang(i,1)=rad2deg(angle(currentline1phn));
        if max([current.a_mang(i) current.b_mang(i) current.c_mang(i)])>Line(i,5)
            cur_viol=cur_viol+1;
        end
    end

    for i=length(Linename)+1 : length(Linename)+length(Cablename)
        DSSCircuit.SetActiveElement(['line.',Cablename{i-length(Linename)}]);  
        curline1=DSSCircuit.ActiveElement.Currents;
        currentline1pha=curline1(1)+1i*curline1(2);
        currentline1phb=curline1(3)+1i*curline1(4);
        currentline1phc=curline1(5)+1i*curline1(6);
        currentline1phn=curline1(7)+1i*curline1(8);
        current.a_mang(i,1)=abs(currentline1pha);
        current.a_ang(i,1)=rad2deg(angle(currentline1pha));
        current.b_mang(i,1)=abs(currentline1phb);
        current.b_ang(i,1)=rad2deg(angle(currentline1phb));
        current.c_mang(i,1)=abs(currentline1phc);
        current.c_ang(i,1)=rad2deg(angle(currentline1phc));
        current.n_mang(i,1)=abs(currentline1phn);
        current.n_ang(i,1)=rad2deg(angle(currentline1phn));
        if max([current.a_mang(i) current.b_mang(i) current.c_mang(i)])>Cable(i-length(Linename),5)
            cur_viol=cur_viol+1;
        end
    end

    DSSCircuit.SetActiveElement('Transformer.Transf2ph1');
    Transf_Power_temp(1,:)=DSSCircuit.ActiveElement.Powers;
    DSSCircuit.SetActiveElement('Transformer.Transf2ph2');
    Transf_Power_temp(2,:)=DSSCircuit.ActiveElement.Powers;
    DSSCircuit.SetActiveElement('Transformer.Transf2ph3');
    Transf_Power_temp(3,:)=DSSCircuit.ActiveElement.Powers;
    Power_transit=sum(Transf_Power_temp(:,7:8));
    S_transit=sqrt(Power_transit(1)^2+Power_transit(2)^2);

    if S_transit>Transf(5)
        cur_viol=cur_viol+1;
    end

    if cur_viol>0
        thermal=0;
    end

    if max(Voltage_temp.Vmax_pu)>1.1
        volt=0;
    end
end
DSSObj.ClearAll;

if ~cnvrgc_fail
    pwr_pv_init = 0; Ptrf = 0; Qtrf = 0;
    Current.a_mang = zeros(77,1); Current.a_ang = zeros(77,1);
    Current.b_mang = zeros(77,1); Current.b_ang = zeros(77,1);
    Current.c_mang = zeros(77,1); Current.c_ang = zeros(77,1);
    Current.n_mang = zeros(77,1); Current.n_ang = zeros(77,1);
    Voltage.VAN_vec = zeros(1,79); Voltage.VBN_vec = zeros(1,79);
    Voltage.VCN_vec = zeros(1,79); Voltage.V1 = zeros(1,79);
    Voltage.V1_mag_pu = zeros(1,79); Voltage.Vmax = zeros(1,79);
    Voltage.Vmax_pu = zeros(1,79);
else
    Voltage = Voltage_temp;
    Ptrf = Power_transit(1);
    Qtrf = Power_transit(2);
    Current = current;
end
    

switch ctrl_type
    case 1
        save([pwd,filesep,'Results',filesep,'UPF_max_pen_pv_alloc_',num2str(idx)],'Voltage','pwr_pv_init','Ptrf','Qtrf','Current');
    case 2
        save([pwd,filesep,'Results',filesep,'Q_P_max_pen_pv_alloc_',num2str(idx)],'Voltage','pwr_pv_init','Ptrf','Qtrf','Current');
    case 3
        save([pwd,filesep,'Results',filesep,'Q_V_max_pen_pv_alloc_',num2str(idx)],'Voltage','pwr_pv_init','Ptrf','Qtrf','Current');

end

end

