clc;
clear all;
close all;

define_constants;

mpc_case = loadcase('case141.m');
base_case = mpc_case;
mpc_case_results = runpf(mpc_case);
number_of_scenarios = 50000;
V = zeros(length(mpc_case.bus),number_of_scenarios);
PD_case = zeros(length(mpc_case.bus),number_of_scenarios);
QD_case = zeros(length(mpc_case.bus),number_of_scenarios);
V(:,1) = mpc_case_results.bus(:,VM);
PD_case(:,1) = base_case.bus(:,PD);
QD_case(:,1) = base_case.bus(:,QD);
i = 1;
while i<=number_of_scenarios
    mult = ones(1,length(mpc_case.bus));
%     Randomly choose the number of loads to be changed
    number_location = randi([1, length(mpc_case.bus)],1,1);
%     Randomly choose the locations of those loads
    locations = randperm(length(mpc_case.bus),number_location);
%     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     matrix
    mult(1,locations) = (1.5+0.5)*rand(length(locations),1)-0.5;
%     Change real and reactive power 
    
    mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ; 
    mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ; 
    [mpc_case_results,success] = runpf(mpc_case);
    if success == 1
        i = i+1;
        PD_case(:,i) = base_case.bus(:,PD).*mult' ;
        QD_case(:,i) = base_case.bus(:,QD).*mult' ;
        V(:,i) = mpc_case_results.bus(:,VM);
    end
end





