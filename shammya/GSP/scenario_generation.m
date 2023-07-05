clc;
clear all;
close all;

define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
number_of_scenario = 8000;
PD_case = zeros(length(mpc_case.bus),number_of_scenario);
QD_case = zeros(length(mpc_case.bus),number_of_scenario);
PD_case(:,1) = base_case.bus(:,PD);
QD_case(:,1) = base_case.bus(:,QD);
avg_pd = mean(base_case.bus(:,PD));
avg_qd = mean(base_case.bus(:,QD));
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);
F= full(sparse(1:nb-1,f,1,nb-1,nb));
T= full(sparse(1:nb-1,t,1,nb-1,nb));
scenario = 1;
nonzero_demands = find(mpc_case.bus(:,PD)>0);
number_location = length(nonzero_demands);
% while scenario <=number_of_scenario
%     mult = ones(1,length(mpc_case.bus));
%     
%     number_location = randi([1, length(mpc_case.bus)],1,1);
%     %     Randomly choose the locations of those loads
%     locations = randperm(length(mpc_case.bus),number_location);
%     %     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     %     matrix
%     mult(1,locations) = (1.5+1.5)*rand(length(locations),1)-1.5;
%     if scenario == 1
%         mult(1,locations) = 1;
%     end
%     mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ;
%     mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ;
%     for l = 1:length(locations)
%         if mpc_case.bus(locations(l),PD) == 0
%             mpc_case.bus(locations(l),PD) = avg_pd*mult(locations(l));
%         end
%         if mpc_case.bus(locations(l),QD) == 0
%             mpc_case.bus(locations(l),QD) = avg_qd*mult(locations(l));
%         end
%         
%     end
%     fprintf('Running Power Flow for Scenario:%d\n',scenario);
%     [mpc_case_results,success] = runpf(mpc_case,opt);
%     if (success == 1)
%         scenario = scenario + 1;
%         PD_case(:,scenario) = mpc_case.bus(:,PD);
%         QD_case(:,scenario) = mpc_case.bus(:,QD);
%     end
%     
% end
while scenario <=number_of_scenario
    mult = ones(1,length(mpc_case.bus));
    mult(1,nonzero_demands) = (1.5+1.5)*rand(length(nonzero_demands),1)-1.5;
    if scenario == 1
        mult(1,nonzero_demands) = 1;
    end
    mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ;
    mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ;
    fprintf('Running Power Flow for Scenario:%d\n',scenario);
    [mpc_case_results,success] = runpf(mpc_case,opt);
    if (success == 1)
        scenario = scenario + 1;
        PD_case(:,scenario) = mpc_case.bus(:,PD);
        QD_case(:,scenario) = mpc_case.bus(:,QD);
    end
end

%%
save('85bus_scenario_only_load.mat','PD_case','QD_case');


% case22_scenario = load('22bus_scenario.mat');