% This code intends to merge the sensor placement algorithm with the
% original optimization problem to see how the sensor placement can still
% generate data for state estimation results. This code also should include
% testing the results against 8000 synthetic cases to generate the
% histogram
clc;
clear all;
close all;

rng(1)
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;


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
cvx_solver SDPT3
cvx_precision high

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
SENSOR_PER_CLUSTER = ceil(nb*0.3);
% U_k = U_Y_new (:,1:number_of_eigenvectors);
U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);

Yuk = Y_DD * U_k;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
Yik=(T'*Y_t+F'*Y_f)*U_Y_new;
for i = 1:ncluster
    nodes = find(Q_L_sym_idx == i);
    [~,sl,al]= greedy_placement(U_k,SENSOR_PER_CLUSTER,[],nodes',G,f,t,4,Yuk,[]);
    sensor_locations = [sensor_locations sl'];
    %     highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end

h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)

%%
% [~,sensor_locations,~]= greedy_placement(U_k,16,[],1:nb,G,f,t,4,Yuk(:,1:number_of_eigenvectors),Yik(:,1:number_of_eigenvectors));
% sensor_locations = [3     5     8     9    12    16    17    19    22    24    25    27    32    34    48    54 58 64 68];

sensor_locations = unique(sensor_locations);
sensor_locations = 1:nb;
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end

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

weights_s = zeros(length(sensor_locations),1);
weights_c = zeros(length(sensor_locations),1);
for i = 1:length(sensor_locations)
    weights_s(i) = norm(Y_DD(i,i));
    weights_c(i) = (norm(Y_DD(i,i))*norm(Y_DD(i,i)));
%     weights_s(i) = vm(i)/abs(sm(i));
%     weights_c(i) = vm(i)/abs(cm(i));
end

% highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
%%
close all 
cvx_clear
cvx_solver Mosek_3
cvx_begin quiet
cvx_precision low
variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable s_est(nb) complex
variable v_est(nb)
variable a_est(nb)
% minimize (65*norm(sm - M*s_est, 2) + 1*norm(vm.*vm-M*(v_est),2) + 4*norm(cm.*cm-M*a_est,2))
minimize (norm(sm - M*s_est, 2) + norm(vm.*vm-M*(v_est),2) + norm(cm.*cm-M*a_est,2))
% minimize (norm(sm - M*s_est, 2))
% minimize norm(min(weights_s)*(sm - M*s_est)) + norm(vm.*vm-M*(v_est),2) + norm(min(weights_c)*(cm.*cm-M*a_est))

subject to
%     tildeW == tildeW';
    s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
    v_est == diag(U_k*tildeW*(U_k)');
%     sum(real(s_est)) <=0.05;
%     sum(imag(s_est)) <=0.05;
    a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
    U_k(1,:)*(tildeW*U_k(1,:)') == 1;
cvx_end


[voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
rmse_angle_single  = rms(angle(vcomplex)-voltage_angle) ;
rmse_voltage_single = rms(abs(vcomplex) - voltage_mag);
fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,rmse_angle_single,rmse_voltage_single);

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

%%
% cvx_solver SDPT3
% all_results = 100*ones(100*100,4);
% c1 = 1;
% scenario = 1;
% for c2 = 1:100
%     for c3 = 1:100
%     cvx_clear
%     
%     cvx_begin quiet
%     cvx_precision high
%     variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
%     variable s_est(nb) complex
%     variable v_est(nb)
%     variable a_est(nb)
%     minimize (c1*norm(sm - M*s_est, 2) + c2/100*norm(vm.*vm-M*(v_est),2) + c3/100*norm(cm.*cm-M*a_est,2))
% 
%     subject to
%         tildeW == tildeW';
%         s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%         v_est == diag(U_k*tildeW*(U_k)');
%         sum(real(s_est)) <=0.05;
%         sum(imag(s_est)) <=0.05;
%         a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
%         U_k(1,:)*(tildeW*U_k(1,:)') == 1;
%     cvx_end
%     [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
%     rmse_angle  = rms(angle(vcomplex)-voltage_angle) ;
%     rmse_voltage = rms(abs(vcomplex) - voltage_mag);
%     
%     if strcmpi('Solved',cvx_status)
%         all_results (scenario,:) = [c2 c3 rmse_angle rmse_voltage];
%         fprintf('Running Case for c2:%d c3:%d\n',c2,c3);
%     end
%     if strcmpi('Inaccurate/Solved',cvx_status)
%         all_results (scenario,:) = [c2 c3 rmse_angle rmse_voltage];
%         fprintf('Running Case for c2:%d c3:%d\n',c2,c3);
%     end
%     scenario = scenario+1;
%     end
% end

%% 
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
