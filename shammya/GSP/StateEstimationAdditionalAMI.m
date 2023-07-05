% clc;
clear all;
close all;

rng(1)
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 0;
bmult=0.25;

casestr = 'case85';
opt=mpoption('out.all',0,'verbose',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
number_of_eigenvectors = 24 ;
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= digraph(f,t,edge_weights);
nonzero_demand = find(mpc_case.bus(:,PD)>0);
ncluster = 1;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.

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
G= graph(f,t,edge_weights);
U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
[U_k,~,~ ] = gram_schemidt(U_k_gft);
Yuk = Y_DD * U_k;

%%
% G= digraph(f,t,edge_weights);
mae_voltage_gsp = [];
mae_angle_gsp = [];
for bmult = 0.0:0.2:1
ami_buses = find(mpc_case.bus(:,PD)>0);
nonami_buses = setdiff(1:nb,ami_buses);
SENSOR_PER_CLUSTER = ceil(length(nonami_buses)*bmult)+length(ami_buses);
[~,sensor_locations,~]= greedy_placement(U_k,SENSOR_PER_CLUSTER,ami_buses,1:nb,G,f,t,4,Yuk,[]);
% number_of_eigenvectors = 16 ;
% sensor_locations = ami_buses;

%%
% sensor_locations = 1:nb;
% sensor_locations = find(mpc_case.bus(:,PD)>0);
% sensor_locations = [1;sensor_locations];
% sensor_locations = 1:nb;
sensor_locations = [1;sensor_locations];
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);

% Choosing the first selected eigenvectors

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

% weights_s = zeros(length(sensor_locations),1);
% weights_c = zeros(length(sensor_locations),1);
% for i = 1:length(sensor_locations)
%     weights_s(i) = norm(Y_DD(i,i));
%     weights_c(i) = (norm(Y_DD(i,i))*norm(Y_DD(i,i)));
% %     weights_s(i) = vm(i)/abs(sm(i));
% %     weights_c(i) = vm(i)/abs(cm(i));
% end



%%
% number_of_eigenvectors = 30 ;
% U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
% [U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);
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

xm=[vm.*vm;conj(sm);cm.*cm];
H_M = [HV;HS;HC];
U_k_kron = kron(conj(U_k),U_k);
H_tilde = H_M*U_k_kron;
Hs_tilde = HS*U_k_kron;

%%
cvx_clear
cvx_solver Mosek
cvx_begin quiet
cvx_precision high
    variable tildeW1(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    minimize norm((xm-H_tilde*vec(tildeW1)))
    subject to
        H_tilde(1,:)*vec(tildeW1) == 1;
cvx_end

if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
    [voltage_mag_gsp,voltage_angle_gsp] = calculate_voltage_phasor(tildeW1,U_k,G);
    mae_angle_gsp  = [mae_angle_gsp mae(angle(vcomplex)-voltage_angle_gsp)] ;
    mae_voltage_gsp = [mae_voltage_gsp mae(abs(vcomplex) - voltage_mag_gsp)];
    fprintf('optval:%f MAE Angle:%f, MAE Voltage: %f\n\n',cvx_optval,mae_angle_gsp,mae_voltage_gsp);
    figure(1)
    subplot(211)
    plot(voltage_mag_gsp)
    hold on
    plot(abs(vcomplex))
    legend('Reconstructed','Original')
    subplot(212)
    plot(voltage_angle_gsp)
    hold on
    plot(angle(vcomplex))
    legend('Reconstructed','Original')
end
end

figure(1)
plot(mae_voltage_gsp,'linewidth',2)
hold on 
plot(mae_angle_gsp,'linewidth',2)
legend('voltage','angle')
%%
% cvx_clear
% cvx_solver Mosek_3
% cvx_precision low
% cvx_begin quiet
%     variable W(nb,nb) complex semidefinite
%     minimize norm((xm - H_M*vec(W))) 
%     subject to
%         H_M(1,:)*vec(W) == 1;
%    
% cvx_end
% 
% if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
%     [voltage_mag_convex,voltage_angle_convex] = calculate_voltage_phasor(W,[],G);
%     mae_angle_convex  = mae(angle(vcomplex)-voltage_angle_convex) ;
%     mae_voltage_convex = mae(abs(vcomplex) - voltage_mag_convex);
%     fprintf('optval:%f MAE Angle:%f, MAE Voltage: %f\n',cvx_optval,mae_angle_convex,mae_voltage_convex);
%     figure(2)
%     subplot(211)
%     plot(voltage_mag_convex)
%     hold on
%     plot(abs(vcomplex))
%     legend('Reconstructed','Original')
%     subplot(212)
%     plot(voltage_angle_convex)
%     hold on
%     plot(angle(vcomplex))
%     legend('Reconstructed','Original')
% end

