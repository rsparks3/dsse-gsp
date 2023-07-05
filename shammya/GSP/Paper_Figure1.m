% clc;
clear all;
close all;

rng(1)
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
bmult=1;

casestr = 'case85';
opt=mpoption('out.all',0,'verbose',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
% bus_number = 200;
% [n,e] = single_feeder_gen(bus_number,bus_number/20,0);
% 
% mpc_case = matpower_fmt(n,e);
% mpc_case = parallel_branch_join(mpc_case);
% mpc_case.branch(:,BR_B) = 0;
% mpc_case.bus(:,BS) = 0;
% mpc_case.branch(1,TAP) = 1;


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

%%

sensor_locations = 1:nb;
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

% Choosing the first selected eigenvectors


number_of_eigenvectors = 50 ;
U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);
id = eye(nb,nb);
HV = [];

for i = 1:nb
    HV=[HV;transpose(kron(id(:,i),id(:,i)))];
end
HV= HV(sensor_locations',:);

HS = [];
for i = 1:nb
    HS=[HS;kron(transpose(id(:,i)),Y_DD(i,:))];
end
HS= HS(sensor_locations',:);

HC = [];
for i = 1:nb
    HC=[HC;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
end
HC= HC(sensor_locations',:);

%%
number_of_scenarios = 50;
voltage_mag_estimated = zeros(nb,number_of_scenarios);
voltage_angle_estimated = zeros(nb,number_of_scenarios);
voltage_mag_mae = zeros(1,number_of_scenarios);
voltage_angle_mae = zeros(1,number_of_scenarios);
for i = 1: number_of_scenarios
    s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
    vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
    s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
    % Add complex noise
    v_noise = ADD_NOISE*random_distribution(nb);
    v = vcomplex+v_noise+1i*v_noise;

    % c = (F'*Y_f-T'*Y_t)./degree_vector*vcomplex;
    CM = Y_bus;
    c = CM*vcomplex;
    c_noise = ADD_NOISE*random_distribution(nb);
    c = c + c_noise + 1i*c_noise;
    cm = abs(c(sensor_locations));

    % Select measurements according to sensor locations
    sm= s(sensor_locations);
    vm = abs(v(sensor_locations));

    xm=[vm.*vm;conj(sm);cm.*cm];
    H_M = [HV;HS;HC];
    U_k_kron = kron(conj(U_k),U_k);
    H_tilde = H_M*U_k_kron;
    Hs_tilde = HS*U_k_kron;

    %%

    cvx_clear
    cvx_solver Mosek_3
    cvx_begin quiet
    cvx_precision high
        variable tildeW1(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
        minimize norm((xm-H_tilde*vec(tildeW1)))
        subject to
            H_tilde(1,:)*vec(tildeW1) == 1;
    cvx_end
    [voltage_mag_new,voltage_angle_new] = calculate_voltage_phasor(tildeW1,U_k,G);
    voltage_mag_estimated(:,i) = voltage_mag_new;
    voltage_angle_estimated(:,i) = voltage_angle_new;
    mae_angle_single_new  = mae(angle(vcomplex)-voltage_angle_new);
    mae_voltage_single_new = mae(abs(vcomplex) - voltage_mag_new);
    voltage_mag_mae (1,i) = mae_voltage_single_new;
    voltage_angle_mae (1,i) = mae_angle_single_new;
    if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
        fprintf('Scenario:%d::optval:%f RMS Angle:%f, RMS Voltage: %f\n',i,cvx_optval,mae_angle_single_new,mae_voltage_single_new);
%         figure(1)
%         subplot(211)
%         plot(voltage_mag_new)
%         hold on
%         plot(abs(vcomplex))
%         legend('Reconstructed','Original')
%         subplot(212)
%         plot(voltage_angle_new)
%         hold on
%         plot(angle(vcomplex))
%         legend('Reconstructed','Original')
    end
end
csvwrite('Paper_Figure1.csv',[transpose(1:nb) mean(voltage_mag_estimated,2) mean(voltage_angle_estimated,2) abs(vcomplex) angle(vcomplex)])
%%
% tic
% cvx_clear
% cvx_solver Mosek_3
% cvx_precision high
% cvx_begin quiet
%     variable W(nb,nb) complex semidefinite
%     minimize norm((xm - H_M*vec(W))) 
%     subject to
%         H_M(1,:)*vec(W) == 1;
%      
% cvx_end
% toc
% if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
%     [voltage_mag2,voltage_angle2] = calculate_voltage_phasor(W,[],G);
%     rmse_angle_single_new  = mae(angle(vcomplex)-voltage_angle2) ;
%     rmse_voltage_single_new = mae(abs(vcomplex) - voltage_mag2);
% 
%     fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle_single_new,rmse_voltage_single_new);
%     figure(2)
%     subplot(211)
%     plot(voltage_mag2)
%     hold on
%     plot(abs(vcomplex))
%     legend('Reconstructed','Original')
%     subplot(212)
%     plot(voltage_angle2)
%     hold on
%     plot(angle(vcomplex))
%     legend('Reconstructed','Original')
% end

% 