% This code generates the plots to show the low rank approximation of the
% paper for mutiple distribution networks, change the distribtion network
% defined by the variable `casestr` and run the same code which will
% generate output in the CSV file
clc;
clear all;
% close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case85.m';
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



[Y_bus,~,~] = makeYbus(mpc_case);
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);
Y_DD = full(Y_bus); % diagonal matrix 

[U_Y, Lambda_Y] = eig(Y_DD,'vector');%./mpcint.baseMVA);
[val_Y, ind_Y]=sort((Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));


[U_Y_new ,~,~]= gram_schemidt(U_Y_new);
%% 
nb = length(mpc_case.bus);
k=nb;
v_k_v = transpose(U_V_new(:,1:k))*Vcomplex(:,1:number_of_scenarios);
v_k_y = transpose(U_Y_new(:,1:k))*Vcomplex(:,1:number_of_scenarios);
diff = abs(v_k_y-v_k_v);
% for i =1:10
%     diff(:,i) = diff(:,i)/val_V(i,1);
% end
fprintf('Min Diff:%f...Max Diff:%f\n', min(abs((v_k_y)-(v_k_v)),[],'all'), max(abs((v_k_y)-(v_k_v)),[],'all'));


%% 

% close all 
% figure
% subplot(311)
% semilogy(mean(abs(v_k_y),2),'linewidth',1.5)
% xlim([0 nb])
% legend('Average GFT Voltage Phasol.,nb99r')
% subplot(312)
% semilogy(std(abs(v_k_y),0,2),'linewidth',1.5)
% xlim([0 nb])
% legend('STD of GFT Voltage Phasor')
% subplot(313)
% semilogy(1./abs(val_Y),'linewidth',1.5)
% legend('1/EigenValues')
% xlim([0 nb])
% 
% 
% close all 
% figure
% semilogy(mean(abs(v_k_y),2),'linewidth',1.5)
% hold on
% semilogy(std(abs(v_k_y),0,2),'linewidth',1.5)
% hold on
% semilogy(1./abs(val_Y),'linewidth',1.5)
% legend('Average GFT Voltage Phasor','STD of GFT Voltage Phasor','1/EigenValues')
% xlim([1 nb])
% % 
% % csvwrite('low_pass.csv',[transpose(1:nb) mean(abs(v_k_v),2) std(abs(v_k_v),0,2) 1./abs(val_Y) ])
% % 
% % 
% close all 
x_axis_data = abs(val_Y)/max(abs(val_Y));
y_axis_mean = mean(abs(v_k_y),2);
y_axis_var = std(abs(v_k_y),0,2);

figure(3)
loglog(x_axis_data,y_axis_mean,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean+y_axis_var,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean-y_axis_var,'linewidth',1.5)
legend('mean','+variance','-variance')
% csvwrite('low_pass_33.csv',[x_axis_data y_axis_mean y_axis_mean+y_axis_var y_axis_mean-y_axis_var ])


% eta_k = zeros(nb-1,1);
% 
% for k = 1:nb-1
% 
%     eta_k(k) = abs(val_Y(k))/abs(val_Y(k+1));
% end
% % eta_k(nb) = val_Y(nb)/min(val_Y(1:nb-1));
% figure()
% subplot(211)
% 
% grid on
% plot((abs(eta_k)),'linewidth',1.5)
% ylabel('eta k')
% xlim([1 nb-1])
% subplot(212)
% 
% grid on
% plot(log(abs(eta_k)),'linewidth',1.5)
% xlim([1 nb-1])
% ylabel('eta k (log scale)')







