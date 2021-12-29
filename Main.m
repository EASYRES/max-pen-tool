clear all
clc

load([pwd,filesep,'PV_allocation',filesep,'Allocation pattern'],'alloc_pv','num_pv','num_allct')

h = waitbar(0,'Please wait...');

ctrl_num = 3 ;   % Number of control types 
steps = ctrl_num*num_allct;

for ctrl_type = 1 : 3
    % 1 for UPF method, 
    % 2 for Q(P) method
    % 3 for Q(V) method

    % Initial overall installed power of PV units (in kw)
    pwr_pv_init = 100; 

    % Uniform distribution of rated power among PV units
    dstr_fctr = 1 / num_pv;

    for i = 1 : num_allct

%       data_write( alloc_pv(i,:), num_pv, i ); Disabled (.dss files have
%       been already created)

        %     Perform power flow analysis to determine the time instant with
        %     maximum voltage
        time_instant = power_flow( num_pv, pwr_pv_init, dstr_fctr, i );

        max_pen_calc( alloc_pv(i,:), num_pv, pwr_pv_init, dstr_fctr, i, ctrl_type, time_instant );
        
        waitbar(((ctrl_type-1)*num_allct+i)/steps)

    end
end

close(h)