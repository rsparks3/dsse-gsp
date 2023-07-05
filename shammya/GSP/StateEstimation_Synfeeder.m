clc;
clear all;
close all;
% rng(1)
define_constants;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.
bus_number = 400;
ADD_NOISE = 0 ;
[n,e] = single_feeder_gen(bus_number);

opt = mpoption('verbose',0,'out.all',0);
mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);
% F= full(sparse(1:nb-1,f,1,nb-1,nb));
% T= full(sparse(1:nb-1,t,1,nb-1,nb));
cvx_solver mosek

%%
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

number_of_eigenvectors = 200;
U_k = U_Y_new (:,1:number_of_eigenvectors);

%% Sensor Placement
ncluster=10;
sensor_per_cluster = 15;
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
Q_sym = full(laplacian_matrix_sym);
[Q_L_sym_whole, ~]=eigs(Q_sym, ncluster,'smallestabs');
Q_L_sym = Q_L_sym_whole;
ammeter_locations = [];
sensor_locations = [];
for i = 1:nb
    Q_L_sym(i,:) =  Q_L_sym(i,:) / sqrt(sum(Q_L_sym(i,:).*Q_L_sym(i,:)));
end
[Q_L_sym_idx, ~] = kmeans(Q_L_sym,ncluster);

Adjacency_clusters_L_sym=zeros(length(Q_L_sym_idx));
for i=1:length(Q_L_sym_idx)
    edges = find(Q_L_sym_idx==Q_L_sym_idx(i)); %from the same group
    Adjacency_clusters_L_sym(i,edges')=1;
end
G= digraph(f,t,edge_weights);
figure(7)
subplot(121)
h = feeder_plot(n,e);
title('Normalized Laplacian and Weighted adjacency matrix')
Yuk = Y_bus * U_Y_new;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
Yik=(T'*Y_t+F'*Y_f)*U_Y_new;
n_eigs = number_of_eigenvectors;
for i = 1:ncluster
    nodes = find(Q_L_sym_idx == i);
    [~,sl,al]= greedy_placement(U_Y_new(:,1:n_eigs),sensor_per_cluster,[],nodes',G,f,t,3,Yuk(:,1:n_eigs),Yik(:,1:n_eigs));
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
% sensor_locations = 1:nb;
sensor_locations = sensor_locations';
sensor_locations = unique(sensor_locations);
number_of_sensors = length(sensor_locations);
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end
subplot(122)
h = feeder_plot(n,e);
title(strcat('Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);


%% State Estimation Routine
[mpc_case_results,success] = runpf(mpc_case,opt);
Uk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
YUk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
UkH_YkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
A=U_k;
C=U_k';
for k = 1:nb
   Uk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end
A=Y_DD*U_k;
for k = 1:nb
   YUk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end

C = U_k'*Y_DD';
for k = 1:nb
   UkH_YkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end



s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
cm = abs(c(sensor_locations));
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3);
slack_bus_voltage = mpc_case.bus(slack_bus,VM);
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
disp('Preprocessing Completed')
%%
% Select measurements according to sensor locations

cvx_clear
cvx_solver SDPT3
cvx_begin
variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable s_est(nb) complex 
variable v_est(nb)
variable a_est(nb)
minimize (15*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

subject to
    tildeW == tildeW';
    for k = 1:nb
        v_est(k) == Uk_UkH(k,:)*vec(transpose(tildeW)); 
        s_est(k) == conj( YUk_UkH(k,:)*vec(transpose(tildeW)) ); 
        a_est(k) == UkH_YkH(k,:)*vec(transpose(tildeW)); 
    end
%     s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%     v_est == diag(U_k*tildeW*(U_k)');
%     sum(real(s_est)) <=0.1;
%     sum(imag(s_est)) <=0.1;
%     a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
    U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_voltage;
cvx_end
fprintf('CSV_STATUS:%s\n',cvx_status);
if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
    [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
    max_angle  = max(abs(angle(vcomplex) - voltage_angle));
    max_voltage  = max(abs(abs(vcomplex) - voltage_mag));
    rmse_angle = rms(angle(vcomplex)-voltage_angle) ;
    rmse_voltage = rms(abs(vcomplex) - voltage_mag);
    fprintf('optval:%f MAX Angle:%f, MAX Voltage: %f, RMS Angle:%f, RMS Voltage: %f\n',cvx_optval,max_angle,max_voltage, rmse_angle,rmse_voltage);
    figure(1)
    subplot(211)
    plot(abs(vcomplex),'linewidth',1.5)
    hold on 
    plot(voltage_mag,'linewidth',1.5)
    legend('Original','Constructed')

    subplot(212)
    plot(angle(vcomplex),'linewidth',1.5)
    hold on 
    plot(voltage_angle,'linewidth',1.5)
    legend('Original','Constructed')
end



%%
tic
cvx_clear
cvx_solver SDPT3
cvx_precision high
cvx_begin quiet
    variable tildeWf(nb,nb) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (15*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

    subject to
        tildeWf == tildeWf';
        s_est == conj(diag(Y_DD*tildeWf));
        v_est == diag(tildeWf);
        a_est == diag(Y_DD*tildeWf*Y_DD');
        tildeWf(1,1) == 1;
cvx_end
fprintf('CSV_STATUS:%s\n',cvx_status);
if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
    toc
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


