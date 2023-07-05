clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case69.m';
mpc_case = loadcase(casestr);
nb = length(mpc_case.bus);

[Y_bus,~,~] = makeYbus(mpc_case);
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);
Y_DD = full(Y_bus)+INCLUDE_LOAD*diag(Y_L); % diagonal matrix 

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));



%%


% number_of_eigenvectors = ceil(1*nb);
  % Get the number of buses
 % Initialize the selection matrix with zeros
[sensor_locations,~] = sensor_placement_matpower(mpc_case);
% [~,sensor_locations] = sensor_placement_matpower(mpc_case);
sensor_locations = [sensor_locations 1];
% Exclude the slack bus measurement
% sensor_locations = sensor_locations(sensor_locations~=1);
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

mpc_case = loadcase(casestr); % load case to get back to the original network
mpc_case_results = runpf(mpc_case); % run the power flow
% Choosing the first selected eigenvectors

number_of_eigenvectors = 8; %chosen arbitarily
U_k = U_Y_new (:,1:number_of_eigenvectors);
s_noise = 0*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;
vm = v(sensor_locations);
v_mag = abs(vm);


f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
% figure
G= graph(f,t);

%% Solve the optimization problem that includes gamma and assumes a constant slack bus
pi_1 = eye(size(nb,nb)) - (U_k*U_k') ;
pi_plus = pi_1(:,1);
pi_minus = pi_1(:,2:end);
M_bar = M(:,2:end);

%% Using the Optimization formulation assuming slack bus voltage = 1 and then gradient descent formula derived from the optimization

cvx_solver mosek
gamma = 0.01:0.01:1;
voltage_error = 10^9*ones(1,length(gamma));
for i = 1:length(gamma)
    cvx_clear
    cvx_begin quiet
        variable v_cvx_1(nb-1) complex
        minimize ((norm(vm-M_bar*v_cvx_1,2)) + gamma(i)*(norm(pi_minus*v_cvx_1+pi_plus,2)))
    cvx_end
voltage_error(i) = norm([1;v_cvx_1]-vcomplex,2);
end

[~,min_gamma_loc] = min(voltage_error);
min_gamma = gamma(min_gamma_loc);
gamma=min_gamma;
gamma = 0.87;
clear v_cvx_1;
cvx_clear
cvx_begin 
    variable v_cvx_2(nb-1) complex
    minimize ((norm(vm-M_bar*v_cvx_2,2)) + gamma*(norm(pi_minus*v_cvx_2+pi_plus,2)))
cvx_end
v_cvx_2 = [1;v_cvx_2];

condition_check = transpose(M_bar)*M_bar + gamma* transpose(pi_minus) * pi_minus ;
rcond(condition_check)
rhs = transpose(M_bar) * vm - gamma* transpose(pi_minus) * pi_plus;
v_bar = [1;pinv(condition_check)* rhs];

figure
subplot(2,1,1)
plot(abs(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((abs(v_cvx_2)),'linewidth',1.5)
hold on 
plot((abs(v_bar)),'linewidth',1.5)
xlim([1 nb])
legend('true','cvx-gamma','gradient-gamma')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Gamma')


subplot(2,1,2)
plot(angle(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((angle(v_cvx_2)),'linewidth',1.5)
hold on 
plot((angle(v_bar)),'linewidth',1.5)
hold off
legend('true','cvx-gamma','gradient-gamma')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Gamma')



%% Using the psueod inverse formulation assuming slack bus voltage = 1
v_pseudo = (U_k(2:end,:))*pinv(M_bar*U_k(2:end,:))*v(sensor_locations);
v_pseudo = [1;v_pseudo];
v_new = (U_k)*pinv(M*U_k)*vm;

figure
subplot(2,1,1)
plot(abs(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((abs(v_new)),'linewidth',1.5)
hold off
xlim([1 nb])
legend('true','pseudo')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Pseudo Inverse')


subplot(2,1,2)
plot(angle(vcomplex),'linewidth',1.5)
hold on 
grid on 

plot((angle(v_new)),'linewidth',1.5)
hold off
legend('true','pseudo')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Pseudo Inverse')


%% Enforce the equality constraint
voltage_error = [];
for gamma = 0.001:0.001:1
    cvx_clear
    cvx_begin quiet
        variable v_new(nb) complex
        minimize (norm(vm-M*v_new,2) + gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new,2))
        subject to 
            v_new(1) == 1;
    cvx_end
voltage_error = [voltage_error norm(v_new-vcomplex,2)];
end
[~,min_gamma_loc] = min(voltage_error);
gamma = 0.001:0.001:1;
min_gamma = gamma(min_gamma_loc);

% min_gamma = 0.87;
cvx_clear
cvx_solver mosek
cvx_begin
    variable v_new_cvx(nb) complex
    minimize (norm(vm-M*v_new_cvx,2) + min_gamma*norm((eye(nb,nb) - (U_k*U_k'))*v_new_cvx,2))
    subject to 
        v_new_cvx(1) == 1;
cvx_end

figure
subplot(2,1,1)
plot(abs(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((abs(v_new_cvx)),'linewidth',1.5)
xlim([1 nb])
legend('true','cvx-gamma-enforce')
xlabel('Bus Numbers')
title('Voltage Magnitude Plot Using Gamma Enforcing Slack bus as constraint')


subplot(2,1,2)
plot(angle(vcomplex),'linewidth',1.5)
hold on 
grid on 
plot((angle(v_new_cvx)),'linewidth',1.5)
hold off
legend('true','cvx-gamma-enforce')
xlim([1 nb])
xlabel('Bus Numbers')
title('Voltage Angle Plot Using Gamma Enforcing Slack bus as constraint')

%%
cvx_clear
cvx_solver mosek
cvx_begin
    variable W(nb,nb) complex semidefinite
    minimize (norm(conj(s(sensor_locations))- diag(M*Y_DD*W*M'),2) + 100*norm(v_mag.*v_mag-diag(M*W*M'),2) + 20*norm(pi_1*W*transpose(pi_1),2))
    subject to 
        W == W';
        W(1,1) == 1;
    
cvx_end

figure 
plot(abs(vcomplex),'linewidth',1.5)
hold on 
plot(abs(sqrt(diag(W))),'linewidth',1.5)
hold off 

%% Enegy in the subspace calculation
base_case = mpc_case;
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



[U_V_l,E_V] = eig(covaraince);
[val_V, ind_V]=sort(diag(E_V),'descend');
U_V_arr = U_V_l(:,ind_V');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_V_arr(:,i))*U_V_arr(:,i));
end
    
U_V_new = U_V_arr*diag(r_vect);

U_tilde = U_V_new(:,1:number_of_eigenvectors);

numerator = norm((eye(nb,nb)-U_k*U_k')*U_tilde,2);
denominator = norm(U_k*U_k'*U_tilde,2);

projection_ratio = numerator^2 / denominator^2;








