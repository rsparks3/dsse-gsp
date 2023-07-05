clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case141';
mpc_case = loadcase(casestr);
base_case = mpc_case;
mpc_case_results = runpf(mpc_case);
number_of_scenarios = 5000;
V = zeros(length(mpc_case.bus),number_of_scenarios);
Vangle = zeros(length(mpc_case.bus),number_of_scenarios);
PD_case = zeros(length(mpc_case.bus),number_of_scenarios);
QD_case = zeros(length(mpc_case.bus),number_of_scenarios);

V(:,1) = mpc_case_results.bus(:,VM);
Vangle(:,1) =  mpc_case_results.bus(:,VA);
PD_case(:,1) = base_case.bus(:,PD);
QD_case(:,1) = base_case.bus(:,QD);
avg_pd = mean(base_case.bus(:,PD));
avg_qd = mean(base_case.bus(:,QD));
i = 1;
while i<number_of_scenarios
    mult = ones(1,length(mpc_case.bus));
%     Randomly choose the number of loads to be changed
    number_location = randi([1, length(mpc_case.bus)],1,1);
%     Randomly choose the locations of those loads
    locations = randperm(length(mpc_case.bus),number_location);
%     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     matrix
    mult(1,locations) = (1.5+1.5)*rand(length(locations),1)-1.5;
%     Change real and reactive power 
    mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ; 
    mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ;
    for l = 1:length(locations)
        if mpc_case.bus(locations(l),PD) == 0
            mpc_case.bus(locations(l),PD) = avg_pd*mult(locations(l));
        end
        if mpc_case.bus(locations(l),QD) == 0
            mpc_case.bus(locations(l),QD) = avg_qd*mult(locations(l));
        end
        
    end
    
    [mpc_case_results,success] = runpf(mpc_case);
    if success == 1
        i = i+1;
        V(:,i) = mpc_case_results.bus(:,VM);
        Vangle(:,i) = mpc_case_results.bus(:,VA);
        PD_case(:,i) = mpc_case.bus(:,PD);
        QD_case(:,i) = mpc_case.bus(:,QD);
    end
end
j=1i;
Vcomplex = zeros(length(mpc_case.bus),number_of_scenarios);
for r = 1:length(mpc_case.bus)
    for c = 1:number_of_scenarios
        Vcomplex (r,c) = V(r,c) * exp(j*Vangle(r,c)/180);
    end
    
end

nb = size(mpc_case.bus,1); %this is the number of buses
nl = size(mpc_case.branch,1); % this is the number of branches

covaraince = Vcomplex * Vcomplex' / number_of_scenarios;


%%
[U_V_l,E_V] = eig(covaraince);
[val_V, ind_V]=sort(diag(E_V),'descend');
U_V_arr = U_V_l(:,ind_V');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_V_arr(:,i))*U_V_arr(:,i));
end
    
U_V_new = U_V_arr*diag(r_vect);

fprintf('Max Difference in estimate for Voltage:%f\n',max(abs(U_V_new * diag(val_V) * transpose(U_V_new) - covaraince),[],'all'));




save(strcat(casestr,'_U_V','.mat'),'U_V_new')








