clc;
clear all;
close all;
rng(1)
define_constants;
SENSOR_PER_CLUSTER = 3;ncluster=8;number_of_eigenvectors=24;
number_of_scenario = 8000;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;

PD_case = load('85bus_scenario_only_load.mat').PD_case;
QD_case = load('85bus_scenario_only_load.mat').QD_case;
nb = length(mpc_case.bus);
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
Q_sym = full(laplacian_matrix_sym);

[U_Q_Y, Lambda_Q_Y] = eig(Q_sym);%./mpcint.baseMVA);
[val_Q_Y, ind_Q_Y]=sort(diag(Lambda_Q_Y),'ascend');

[Q_L_sym_whole, ~]=eigs(Q_sym, ncluster,'smallestabs');
Q_L_sym = Q_L_sym_whole;
ammeter_locations = [];
sensor_locations = [];
for i = 1:nb
    Q_L_sym(i,:) =  Q_L_sym(i,:) / sqrt(sum(Q_L_sym(i,:).*Q_L_sym(i,:)));
end
[Q_L_sym_idx, ~] = kmeans(Q_L_sym,ncluster);
%%
Adjacency_clusters_L_sym=zeros(length(Q_L_sym_idx));
for i=1:length(Q_L_sym_idx)
    edges = find(Q_L_sym_idx==Q_L_sym_idx(i)); %from the same group
    Adjacency_clusters_L_sym(i,edges')=1;
end
G= digraph(f,t,edge_weights);
% figure(7)
% subplot(121)
% h = plot(G);
% title('Normalized Laplacian and Weighted adjacency matrix')
U_k = U_Y_new (:,1:number_of_eigenvectors);
Yuk = Y_bus * U_Y_new;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
Yik=(T'*Y_t+F'*Y_f)*U_Y_new;
for i = 1:ncluster
    nodes = find(Q_L_sym_idx == i);
%     nodes = intersect(nonzero_demand,nodes);
    [~,sl,al]= greedy_placement(U_k,SENSOR_PER_CLUSTER,[],nodes',G,f,t,4,Yuk(:,1:number_of_eigenvectors),Yik(:,1:number_of_eigenvectors));
    sensor_locations = [sensor_locations sl'];
    %     highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


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

%%

actual_angle = zeros(nb,number_of_scenario);
actual_voltage = zeros(nb,number_of_scenario);

estimated_angle = zeros(nb,number_of_scenario);
estimated_voltage = zeros(nb,number_of_scenario);
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
    cvx_solver SDPT3
    cvx_begin quiet
        variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
        variable s_est(nb) complex
        variable v_est(nb)
        variable a_est(nb)
        minimize (5*norm(sm - M*s_est, 2) + 4*norm(vm.*vm-M*(v_est),2) + 4*norm(cm.*cm-M*a_est,2) )

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
        estimated_angle(:,scenario) = voltage_angle;
        estimated_voltage(:,scenario) = voltage_mag;
        actual_voltage(:,scenario) = abs(vcomplex);
        actual_angle(:,scenario) = angle(vcomplex);
        fprintf('Scenario:%d optval:%f\n',scenario,cvx_optval);
    
    else
        failed_scenario = [failed_scenario scenario];
    end    
end

save(strcat(casestr,'_clustering_scenario.mat'),'failed_scenario','estimated_angle','actual_angle','estimated_voltage','actual_voltage')



% save('85Bus_85Sensors.mat','rmse_angle','rmse_voltage','max_angle','max_voltage','failed_scenario')







