clc;
clear all;
close all;
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
nb = length(mpc_case.bus);
k=nb;
v_k_v = transpose(U_V_new(:,1:k))*Vcomplex(:,1:20);
v_k_y = transpose(U_Y_new(:,1:k))*Vcomplex(:,1:20);
diff = abs(v_k_y-v_k_v);
% for i =1:10
%     diff(:,i) = diff(:,i)/val_V(i,1);
% end
fprintf('Min Diff:%f...Max Diff:%f\n', min(abs((v_k_y)-(v_k_v)),[],'all'), max(abs((v_k_y)-(v_k_v)),[],'all'));


%% 

figure

for i = 1:20
  semilogy(abs(v_k_y(:,i)))
  hold on
end
hold off


%%
number_of_sensors = ceil(1*nb-30);
number_of_eigenvectors = ceil(nb)-5;
% number_of_eigenvectors = ceil(1*nb);
  % Get the number of buses
 % Initialize the selection matrix with zeros
sensor_locations = randperm(nb,number_of_sensors);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

mpc_case = loadcase(casestr); % load case to get back to the original network
mpc_case_results = runpf(mpc_case); % run the power flow
% Choosing the first selected eigenvectors
U_k = U_Y_new (:,1:number_of_eigenvectors);
s_noise = 0*random_distribution(nb); % A customized function to generate random values 
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
% Add complex noise
v_noise = 0*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

%%
% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = vcomplex(sensor_locations);

cvx_solver mosek
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
     minimize (norm(sm - M*conj(diag(Y_DD*U_k*tildeW*U_k')), 2) + norm(vm.*vm-diag(M*U_k*tildeW*(U_k)'*M'),2) + 0*trace(tildeW * abs(diag(val_Y(1:number_of_eigenvectors)))))  
    
    subject to
        tildeW == tildeW';

cvx_end


W = U_k * tildeW * U_k';

figure 
subplot(211)
plot(abs(vcomplex))
hold on 
grid on 
plot(sqrt(abs(diag(W))))
hold on 
% plot(abs(v_recon_1))
% plot(abs(v_recon_2))
legend('true','Relaxation')
xlabel('Bus Numbers')
title('Reconstruction with Tilde W')
hold off

f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

G= graph(f,t);

voltage_angle= zeros(nb,1);
for i = 2:nb
    path = shortestpath(G,1,i);
    angle_sum = 0;
    for j = 1:length(path)-1
        angle_sum = angle_sum - angle(W(path(j),path(j+1)));
    end
    voltage_angle(i,1) = angle_sum;
    
end


subplot(212)
plot(mpc_case_results.bus(:,VA),'linewidth',1.5)
hold on 
plot(voltage_angle*180/pi,'linewidth',1.5)
grid on
hold off
xlabel('Bus Number')
ylabel('Voltage Angle (degree)')
xlim([1 nb])
xlabel('Bus Number')
legend('Original','Convex Relaxation')
title('Voltage Angle')


%%












