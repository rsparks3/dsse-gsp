% clc;
clear all;
close all;
define_constants;
casestr = 'case69';
number_of_eigenvectors = 40;
ADD_NOISE = 1;
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;

PD_case = base_case.bus(:,PD);
QD_case = base_case.bus(:,QD);

nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);

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
mpc_case_results = runpf(mpc_case,opt); % run the power flow

sensor_locations = 1:nb;
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end


%%
s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
%%
number_of_eigenvectors = 30;
U_k = U_Y_new (:,1:number_of_eigenvectors);
tic
cvx_clear
cvx_solver SDPT3
cvx_precision high
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (35*norm(sm - M*s_est, 2) + 25*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

    subject to
        tildeW == tildeW';
        s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
        v_est == diag(U_k*tildeW*(U_k)');
        a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
        U_k(1,:)*(tildeW*U_k(1,:)') == 1;
        sum(real(s_est)) <=0.1; 
        sum(imag(s_est)) <=0.1; 
cvx_end
toc
fprintf('cvx_status:%s',cvx_status)
if strcmpi('Solved',cvx_status)||strcmpi('Inaccurate/Solved',cvx_status)
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
    rmse_angle  = rms(angle(vcomplex)-voltage_angle) ;
    rmse_voltage = rms(abs(vcomplex) - voltage_mag);
    fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);
    
end
%%
tic
cvx_clear
cvx_solver SDPT3
cvx_precision high
cvx_begin 
    variable tildeWf(nb,nb) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (5*norm(sm - M*s_est, 2) + 4*norm(vm.*vm-M*(v_est),2) + 4*norm(cm.*cm-M*a_est,2) )

    subject to
        tildeWf == tildeWf';
        s_est == conj(diag(Y_DD*tildeWf));
        v_est == diag(tildeWf);
%         a_est == diag(Y_DD*tildeWf*Y_DD');
        tildeWf(1,1) == 1;
cvx_end
toc
fprintf('cvx_status:%s',cvx_status)
if strcmpi('Solved',cvx_status)
    [voltage_mag_a,voltage_angle_a] = calculate_voltage_phasor(tildeWf,U_k,G);
    rmse_angle  = rms(angle(vcomplex)-voltage_angle_a) ;
    rmse_voltage = rms(abs(vcomplex) - voltage_mag_a);
    fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle,rmse_voltage);
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
end








