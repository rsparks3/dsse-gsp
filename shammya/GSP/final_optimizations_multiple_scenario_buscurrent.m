 % clc;
clear all;
% close all;
% rng(1)   
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
number_of_scenario = 1;
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
%% Calcualtion for Y matrix
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_DD = full(Y_bus);

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
mpc_case_results = runpf(mpc_case); % run the power flow
%% 

number_of_eigenvectors = 23;
%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);

[~,sensor_locations] = sensor_placement_matpower(mpc_case);
ammeter = load('ammeter.mat');

sensor_to_exclude=[76 56 43 22 85 77 22 47 11 13 80 69 40 45 38 67 37 23 66 84 82 78 15 11 54 55 51 18 20 71 75 16 39 34 36 42 74 79 72 59 11 62];
for s_e = 1:length(sensor_to_exclude)
    sensor_locations = sensor_locations(sensor_locations~=sensor_to_exclude(s_e));
end
sensor_locations = [sensor_locations [32 12 34 5 68 27 54 9 48 64 16 8 25 58 22 19 22 3]];
% sensor_locations=[61 57 63 29 26 31 10 1 25 73 65 69 83 14 11 41 53 50 21 19 20 46 32 40];
ammeter_locations  = ammeter.ammeter_location;
ammeter_locations = ammeter_locations(ammeter_locations~=45);
ammeter_locations = ammeter_locations(ammeter_locations~=46);
ammeter_locations = ammeter_locations(ammeter_locations~=13);
ammeter_locations = ammeter_locations(ammeter_locations~=61);
ammeter_locations = ammeter_locations(ammeter_locations~=78);
ammeter_locations = ammeter_locations(ammeter_locations~=81);
ammeter_locations = ammeter_locations(ammeter_locations~=76);
ammeter_locations = ammeter_locations(ammeter_locations~=74);
ammeter_locations = ammeter_locations(ammeter_locations~=11);
ammeter_locations = ammeter_locations(ammeter_locations~=31);
ammeter_locations = [ammeter_locations 48 21 43 53 79 22 39];

M_a =zeros(length(ammeter_locations),nb-1);
for r = 1:length(ammeter_locations)
    M_a (r,ammeter_locations(r)) = 1; % Put 1 in sensor locations
end
sensor_locations = sensor_locations';
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

% Choosing the first selected eigenvectors



figure(1)
h = plot(G);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


max_angle_difference_scenario = zeros(1,number_of_scenario);
max_voltage_difference_scenario = zeros(1,number_of_scenario);

for scenario = 1:number_of_scenario
    mult = ones(1,length(mpc_case.bus));

    number_location = randi([1, length(mpc_case.bus)],1,1);
%     Randomly choose the locations of those loads
    locations = randperm(length(mpc_case.bus),number_location);
%     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     matrix
    mult(1,locations) = (1.5+1.5)*rand(length(locations),1)-1.5;
    if scenario == 1
        mult(1,locations) = 1;
    end
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
    
    [mpc_case_results,success] = runpf(mpc_case,opt);
    
    PD_case(:,scenario) = mpc_case.bus(:,PD);
    QD_case(:,scenario) = mpc_case.bus(:,QD);
    s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values 
    vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
    s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
    v_noise = ADD_NOISE*random_distribution(nb);
    v = vcomplex+v_noise+1i*v_noise;

    c = T'*Y_t*vcomplex;
    c_noise = ADD_NOISE*random_distribution(nb);
    c = c + c_noise + 1i*c_noise;
    cm = abs(c(sensor_locations));


    % Select measurements according to sensor locations
    sm= s(sensor_locations);
    vm = abs(vcomplex(sensor_locations));
    cvx_clear
    cvx_solver mosek
    cvx_precision high
    cvx_begin quiet
        cvx_precision high
        variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
        variable s_est(nb) complex
        variable v_est(nb)
        variable a_est(nb) 
        minimize (1*norm(sm - M*s_est, 2) + 30*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )  

        subject to
    %         sm == M*s_est;
            tildeW == tildeW'; 
            s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
            v_est == diag(U_k*tildeW*(U_k)');
            sum(real(s_est)) <=0.1; 
            sum(imag(s_est)) <=0.1; 
            a_est == diag((T'*Y_t)*U_k*tildeW*U_k'*(T'*Y_t)');
            U_k(1,:)*(tildeW*U_k(1,:)') == 1;
    cvx_end


    [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
    max_angle_diff = norm(angle(vcomplex) - voltage_angle,2);
    max_voltage_diff = norm(abs(vcomplex) - voltage_mag,2);
    max_angle_difference_scenario(1,scenario) = max_angle_diff;
    max_voltage_difference_scenario (1,scenario) = max_voltage_diff;
    fprintf('Maximum Angle Diff : %f & Magnitude Difference: %f and optval:%f\n',max_angle_diff,max_voltage_diff,cvx_optval);
end

figure()
subplot(211)
histogram(max_voltage_difference_scenario)
title('Voltage Magnitude Estimation RMSE ')
subplot(212)
histogram(max_angle_difference_scenario)
title('Voltage Angle Estimation RMSE ')



save('8000scenarios_heuristic_to_buscurrent.mat','max_voltage_difference_scenario','max_angle_difference_scenario','PD_case','QD_case')








