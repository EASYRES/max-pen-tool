function [ time_instant ] = power_flow( num_pv, pwr_pv_init, dstr_fctr, idx )

load([pwd,filesep,'Data',filesep,'Consumption and generation data'])

load([pwd,filesep,'Data',filesep,'Network Data'],'Busname','Cablename','Linename','PQLoadname');

for i = 1 : num_pv
    PV_name(i) = {['PV',num2str(i)]};
end

Pmppt=dstr_fctr*pwr_pv_init*Pprofile*ones( 1 , num_pv );

Voltage.V1_mag_pu=[];

for t = 1 : length(PQLoad(:,1))
% Initialization of the OpenDSS
DSSObj = actxserver('OpenDSSEngine.DSS');
Start = DSSObj.Start(0);
DSSCircuit = DSSObj.ActiveCircuit;
DSSText = DSSObj.Text;
DSSObj.ClearAll;
% path of the .dss file
DSSText.Command = ['compile ',pwd,filesep,'Data',filesep,'Examined_Network_',num2str(idx),'.dss'];

% Tap Determination
DSSText.Command = ['Transformer.Transf2ph1.Taps=[',num2str(0.95),',1]'];
DSSText.Command = ['Transformer.Transf2ph2.Taps=[',num2str(0.95),',1]'];
DSSText.Command = ['Transformer.Transf2ph3.Taps=[',num2str(0.95),',1]'];


for i = 1 : length(PQLoadname)
    DSSCircuit.SetActiveElement(['load.',PQLoadname{i}]);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(PQLoad(t,i)));
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(imag(PQLoad(t,i)));
end

for i = 1 : length(PV_name)
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'a']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(t,i))/3);
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);
    
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'b']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(t,i))/3);
    DSSCircuit.ActiveElement.Properties('kvar').Val=num2str(0);
    
    DSSCircuit.SetActiveElement(['generator.',PV_name{i},'c']);  
    DSSCircuit.ActiveElement.Properties('kW').Val=num2str(real(Pmppt(t,i))/3);
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
end
    
Voltage.V1_mag_pu=[Voltage.V1_mag_pu;Voltage_temp.V1_mag_pu];
    
end

temp = max(Voltage.V1_mag_pu');
[~, time_instant] = max(temp);

DSSObj.ClearAll;
end