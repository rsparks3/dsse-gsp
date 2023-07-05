clc
clear
close all;
all_distributions = load('distribution.mat').all_distributions;
all_data = load('all_data.mat').all_data;
data_file_list = load('data_file_list.mat').data_files;
define_constants
mpc = loadcase('case85.m');
P_D = mpc.bus(:,PD);
Q_D = mpc.bus(:,QD);
load_profile_map = randi([1 length(data_file_list)],length(mpc.bus),1);
save('load_profile_map.mat','load_profile_map');
power_flow = cell(1,50000);
s= 1;
while s <= 50000
    for i = 1:length(mpc.bus)
        max_value = max(all_data(data_file_list(load_profile_map(i))));
        mult = random(all_distributions(data_file_list(load_profile_map(i))));
        mult = mult/max_value;
        mpc.bus(i,PD) = mult * P_D(i);
        mpc.bus(i,QD) = mult * Q_D(i);
    end
    [mpc_case_results,success] = runpf(mpc);
    if success == 1
        power_flow{s} = mpc_case_results;
        s = s + 1;
    end
    
end

save('power_flow.mat','power_flow','-v7.3')