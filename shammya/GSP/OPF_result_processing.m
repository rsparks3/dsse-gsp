clc;
clear all;
close all;
define_constants;
opf_cases = [200,400,600,800,2000];

final_result= zeros(length(opf_cases),12);
for i = 1:length(opf_cases)
    bus = opf_cases(i);
    final_result(i,1) = bus;
    load(strcat('OPF_GSP_',string(bus)));
    mpc_case_new=mpc_case;
    % total_kw_demand = sum(mpc_case.bus(:,PD));
    % total_kvar_demand = sum(mpc_case.bus(:,QD));
    total_demand = sum(mpc_case.bus(:,PD)) + 1i*sum(mpc_case.bus(:,QD));
    final_result(i,2) = total_demand;
    total_gen = mpc_case_new.baseMVA*sum(pg) + 1i*mpc_case_new.baseMVA*sum(qg);
    final_result(i,3) = total_gen;
    cost_gsp = sum(Cs.*pg(generation_bus).*pg(generation_bus))*mpc_case_new.baseMVA*mpc_case_new.baseMVA;
    final_result(i,4) = cost_gsp;
    
    slack_bus = find(mpc_case_new.bus(:,BUS_TYPE)==3);

    % mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
    % mpc_case_new.gen = mpc_case_new.gen(1,:);
    mpc_case_new.gen(slack_bus,VG) = voltage_mag(1);
    gen = mpc_case_new.gen;
    all_gen = gen;
    for j = 2:length(generation_bus)
        all_gen = [all_gen;gen];
    end
    all_gen(:,1) = generation_bus;
    all_gen(2:end,PG) = mpc_case_new.baseMVA*pg(generation_bus(2:end));
    all_gen(2:end,QG) = mpc_case_new.baseMVA*qg(generation_bus(2:end));
    all_gen(2:end,QMAX) = mpc_case_new.baseMVA*qg(generation_bus(2:end));
    all_gen(2:end,QMIN) = mpc_case_new.baseMVA*qg(generation_bus(2:end));
    mpc_case_new.gen = all_gen;
    % mpc_case_new.bus(generation_bus(2:end),PD) = mpc_case_new.bus(generation_bus(2:end),PD) - mpc_case_new.baseMVA*pg(generation_bus(2:end));
    % mpc_case_new.bus(generation_bus(2:end),QD) = mpc_case_new.bus(generation_bus(2:end),QD) - mpc_case_new.baseMVA*qg(generation_bus(2:end));


    mpc_case_new_results = runpf(mpc_case_new);
    new_voltage = mpc_case_new_results.bus(:,VM);
    [loss,fchg,tchg] = get_losses(mpc_case_new_results);
    vcomplex = mpc_case_new_results.bus(:,VM).*exp(1j*mpc_case_new_results.bus(:,VA)*pi/180); % Preparing the complex voltages
    [Y_bus,~,~] = makeYbus(mpc_case_new);
    s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
    final_result(i,5) = mean(abs(new_voltage-voltage_mag));
    final_result(i,6) = std(abs(new_voltage-voltage_mag));
    final_result(i,7) = telapsed_gsp;
    
    
    load(strcat('OPF_NOGSP_',string(bus)));
    mpc_case_new=mpc_case;
    final_result(i,8) = mpc_case_new.baseMVA*sum(pg_convex) + 1i*mpc_case_new.baseMVA*sum(qg_convex);
    final_result(i,9) = sum(Cs.*pg_convex(generation_bus).*pg_convex(generation_bus))*mpc_case_new.baseMVA*mpc_case_new.baseMVA;
    mpc_case_new=mpc_case;
    slack_bus = find(mpc_case_new.bus(:,BUS_TYPE)==3);

    % mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
    % mpc_case_new.gen = mpc_case_new.gen(1,:);
    mpc_case_new.gen(slack_bus,VG) = voltage_mag_convex(1);
    gen = mpc_case_new.gen;
    all_gen = gen;
    for j = 2:length(generation_bus)
        all_gen = [all_gen;gen];
    end
    all_gen(:,1) = generation_bus;
    all_gen(2:end,PG) = mpc_case_new.baseMVA*pg_convex(generation_bus(2:end));
    all_gen(2:end,QG) = mpc_case_new.baseMVA*qg_convex(generation_bus(2:end));
    all_gen(2:end,QMAX) = mpc_case_new.baseMVA*qg_convex(generation_bus(2:end));
    all_gen(2:end,QMIN) = mpc_case_new.baseMVA*qg_convex(generation_bus(2:end));
    mpc_case_new.gen = all_gen;
    % mpc_case_new.bus(generation_bus(2:end),PD) = mpc_case_new.bus(generation_bus(2:end),PD) - mpc_case_new.baseMVA*pg(generation_bus(2:end));
    % mpc_case_new.bus(generation_bus(2:end),QD) = mpc_case_new.bus(generation_bus(2:end),QD) - mpc_case_new.baseMVA*qg(generation_bus(2:end));


    mpc_case_new_results = runpf(mpc_case_new);
    new_voltage = mpc_case_new_results.bus(:,VM);
    vcomplex = mpc_case_new_results.bus(:,VM).*exp(1j*mpc_case_new_results.bus(:,VA)*pi/180); % Preparing the complex voltages
    [Y_bus,~,~] = makeYbus(mpc_case_new);
    s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
    final_result(i,10) = mean(abs(new_voltage-voltage_mag_convex));
    final_result(i,11) = std(abs(new_voltage-voltage_mag_convex));
    final_result(i,12) = telapsed;
    

end
csvwrite('C:\Users\sssaha\Dropbox (ASU)\Shammya\Projects\Neptune 2.0\GSP\OPF_results.csv',final_result)
