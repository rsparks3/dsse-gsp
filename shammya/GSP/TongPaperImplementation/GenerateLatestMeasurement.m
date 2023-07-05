define_constants
new_data_file_list= load('new_data_file_list.mat').new_data_files;
load_profile_map = load('load_profile_map.mat').load_profile_map;
all_data_2020 = load('all_data_2020.mat').all_data_2020;
ammeter_location = load('ammeter.mat').ammeter_location;
mpc = loadcase('case85.m');
P_D = mpc.bus(:,PD);
Q_D = mpc.bus(:,QD);
input_data_2020 = zeros(8000,84*2+length(ammeter_location));
output_data_2020 = zeros(8000,85*2);
[~,Yf,~] = makeYbus(loadcase('case85.m'));
for s = 1:8000
    for i = 1:length(mpc.bus)
        demand = all_data_2020(new_data_file_list(load_profile_map(i)));
        
        mult = demand(s)/max(demand);
        if s == 1
            mult = 1;
        end
       
        mpc.bus(i,PD) = mult * P_D(i);
        mpc.bus(i,QD) = mult * Q_D(i);
    end
    [mpc_result,success] = runpf(mpc);
    if success == 1
        pf = mpc_result.branch(:,PF)+random_distribution(84);
        qf = mpc_result.branch(:,QF)+random_distribution(84);
        current = Yf*mpc_result.bus(:,VM);
        current = current(ammeter_location) ;
        current = abs(current) + random_distribution(length(current),0.001);
        input_data_2020(s,:) = [pf' qf' current'];
        vmagnitude = mpc_result.bus(:,VM);
        vangle = mpc_result.bus(:,VA);
        output_data_2020(s,:) = [vmagnitude' vangle'];
    end
end

save('input_data_2020.mat','input_data_2020');
save('output_data_2020.mat','output_data_2020');



