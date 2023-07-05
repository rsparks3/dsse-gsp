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

PD_case = load('85bus_scenario.mat').PD_case;
QD_case = load('85bus_scenario.mat').QD_case;
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);
F= full(sparse(1:nb-1,f,1,nb-1,nb));
T= full(sparse(1:nb-1,t,1,nb-1,nb));
cvx_solver mosek
% cvx_precision 
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

number_of_eigenvectors = 24;
%U_V_new = load(strcat(casestr,'_U_V','.mat')).U_V_new;
U_k = U_Y_new (:,1:number_of_eigenvectors);

% Sensor_locations for 85 bus
% sensor_locations = [44 46 83 19 41 61 53 50 73 65 70 14 29 33 26 7];
sensor_locations = [7    10    14    19    21    26    29    35    41    44    46    58    61    65    73    81];
% 32 Sensors for the 85 bus test case 
% sensor_locations = [44    46    83    19    41    61    53    50    73    65    70    14 58     1    27    11    21    35    18    52    13     8    64    49];
% sensor_locations = 1:85;
% sensor locations for 69 bus using clustering 
% sensor_locations = [34 33 49 1 64 62 40 45 57 56 68 66 14 13 26 24];
% Sensor location for 141 bus clustering 
% sensor_locations=[104 115 39 6 129 124 31 25 58 67 33 5 47 46 81 51];
% Sensor location for 141 bus overall placement 
% sensor_locations = [104    58    33    81   115    67    39   108   123    25    73   129   131    51    21    61];


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


number_of_scenario = 8000;
max_angle = zeros(1,number_of_scenario);
max_voltage = zeros(1,number_of_scenario);

rmse_angle = zeros(1,number_of_scenario);
rmse_voltage = zeros(1,number_of_scenario);
failed_scenario = [];
for scenario = 1:number_of_scenario
    mpc_case.bus(:,PD) = PD_case(:,scenario);
    mpc_case.bus(:,QD) = QD_case(:,scenario);
    [mpc_case_results,success] = runpf(mpc_case,opt);
    
    s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
    vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
    s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
    v_noise = ADD_NOISE*random_distribution(nb);
    v = vcomplex+v_noise+1i*v_noise;
    
    c = Y_DD*vcomplex;
    c_noise = ADD_NOISE*random_distribution(nb);
    c = c + c_noise + 1i*c_noise;
    cm = abs(c(sensor_locations));
    
    
    % Select measurements according to sensor locations
    sm= s(sensor_locations);
    vm = abs(vcomplex(sensor_locations));
    cvx_clear
    cvx_solver mosek
    cvx_begin quiet
        variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
        variable s_est(nb) complex
        variable v_est(nb)
        variable a_est(nb)
        minimize (1*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

        subject to
        %         sm == M*s_est;
        tildeW == tildeW';
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        sum(real(s_est)) <=0.1;
        sum(imag(s_est)) <=0.1;
        a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
    cvx_end
    
    
    if strcmpi('Solved',cvx_status)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        rmse_angle (1,scenario) = rms(angle(vcomplex)-voltage_angle) ;
        rmse_voltage (1,scenario) = rms(abs(vcomplex) - voltage_mag);
        max_angle (1,scenario) = max(abs(angle(vcomplex) - voltage_angle));
        max_voltage (1,scenario) = max(abs(abs(vcomplex) - voltage_mag));
        fprintf('Scenario:%d optval:%f RMS Angle:%f, RMS Voltage: %f\n',scenario,cvx_optval,rmse_angle (1,scenario),rmse_voltage (1,scenario));
    
    else
        failed_scenario = [failed_scenario scenario];
    end    
end

figure()
subplot(211)
histogram(rmse_voltage,'Normalization','probability')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
title('Voltage Magnitude Estimation RMSE','Fontsize',15)
subplot(212)
histogram(rmse_angle,'Normalization','probability')
title('Voltage Angle Estimation RMSE ','Fontsize',15)
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;



% save('85Bus_85Sensors.mat','rmse_angle','rmse_voltage','max_angle','max_voltage','failed_scenario')


%%
scenario = 5229;
mpc_case.bus(:,PD) = PD_case(:,scenario);
mpc_case.bus(:,QD) = QD_case(:,scenario);
[mpc_case_results,success] = runpf(mpc_case,opt);

s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));


% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
cvx_clear
cvx_solver mosek
cvx_begin quiet
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (1*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

    subject to
            sm == M*s_est;
    tildeW == tildeW';
    s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
    v_est == diag(U_k*tildeW*(U_k)');
    sum(real(s_est)) <=0.1;
    sum(imag(s_est)) <=0.1;
    a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
    U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end

[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);


figure
subplot(211)
plot(voltage_mag)
hold on
plot(abs(vcomplex))
legend('Reconstructed','Original')
subplot(212)
plot(voltage_angle)
hold on
plot(angle(vcomplex))
legend('Reconstructed','Original')




