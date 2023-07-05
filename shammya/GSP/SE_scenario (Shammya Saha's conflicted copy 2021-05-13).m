clc;
clear all;
close all;
rng(1)
define_constants;
ncluster=1;
number_of_eigenvectors=24;
number_of_scenario = 8000;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;

PD_case = load('85bus_scenario.mat').PD_case;
QD_case = load('85bus_scenario.mat').QD_case;
nb = length(mpc_case.bus);
SENSOR_PER_CLUSTER = ceil(nb*0.25);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
 
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
G= graph(f,t,edge_weights);
F= full(sparse(1:G.numedges,f,1,G.numedges,G.numnodes));
T= full(sparse(1:G.numedges,t,1,G.numedges,G.numnodes));
degree_matrix = diag(degree(G,1:nb));
degree_vector = diag(degree_matrix);
degree_vector(degree_vector>1) = 2 ;
adjacency_matrix = adjacency(G,'weighted');
% degree_matrix = diag(sum(adjacency_matrix,2));
% laplacian_matrix = degree_matrix-full(adjacency_matrix);
laplacian_matrix_sym = speye(nb,nb) - (degree_matrix)^(-0.5)*adjacency_matrix*(degree_matrix)^(-0.5);
% laplacian_matrix_sym = speye(nb,nb) - (degree_matrix)^(-1)*adjacency_matrix;

%%
G= digraph(f,t,edge_weights);
% figure(7)
% subplot(121)
% h = plot(G);
% title('Normalized Laplacian and Weighted adjacency matrix')
U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
[U_k,~,~ ] = gram_schemidt(U_k_gft);
Yuk = Y_bus * U_k;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
% Yik=(T'*Y_t+F'*Y_f)*U_Y_new;
nodes = 1:nb;
[~,sensor_locations,al]= greedy_placement(U_k,SENSOR_PER_CLUSTER,[],nodes',G,f,t,4,Yuk,[]);

subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)

sensor_locations = [1;sensor_locations];
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

%%

actual_angle = zeros(nb,number_of_scenario);
actual_voltage = zeros(nb,number_of_scenario);

estimated_angle = zeros(nb,number_of_scenario);
estimated_voltage = zeros(nb,number_of_scenario);
failed_scenario = [];

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


H_M = [HV;HS;HC];
U_k_kron = kron(conj(U_k),U_k);
H_tilde = H_M*U_k_kron;
Hs_tilde = HS*U_k_kron;

%%
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
    xm=[vm.*vm;conj(sm);cm.*cm];
    cvx_clear
    cvx_solver Mosek_3
    cvx_begin quiet
    cvx_precision high
        variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
        minimize norm((xm-H_tilde*vec(tildeW)))
        subject to
            H_tilde(1,:)*vec(tildeW) == 1;
    cvx_end
    
    
    if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        estimated_angle(:,scenario) = voltage_angle;
        estimated_voltage(:,scenario) = voltage_mag;
        actual_voltage(:,scenario) = abs(vcomplex);
        actual_angle(:,scenario) = angle(vcomplex);
        fprintf('Scenario:%d optval:%f\n',scenario,cvx_optval);
    
    else
        failed_scenario = [failed_scenario scenario];
    end    
end

save(strcat(casestr,'_25p_mosek.mat'),'failed_scenario','estimated_angle','actual_angle','estimated_voltage','actual_voltage','sensor_locations')



% save('85Bus_85Sensors.mat','rmse_angle','rmse_voltage','max_angle','max_voltage','failed_scenario')







