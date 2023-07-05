% clc;
clear all;
% close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
mpc_case = loadcase(casestr); % load case to get back to the original network
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= digraph(f,t,edge_weights);

ncluster = 8;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.
cvx_solver mosek
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
mpc_case_results = runpf(mpc_case); % run the power flow


%%
G= graph(f,t);
adjacency_matrix = adjacency(G);


%% Un-normalized Graph Laplacian and non-weighted adjacency matrix 
close all 
G= graph(f,t,edge_weights);
degree_matrix = diag(degree(G,1:nb));
adjacency_matrix = adjacency(G);
laplacian_matrix = degree_matrix-full(adjacency_matrix);
ammeter_locations = [];
sensor_locations = [];
Q = full(laplacian_matrix);
[QV, ~]=eigs(Q, ncluster,'smallestabs');
[Qidx, ~] = kmeans(QV(:,1:ncluster),ncluster);
%plot the adjacency:
Adjacency_clusters=zeros(length(Qidx));
for i=1:length(Qidx)
    edges = find(Qidx==Qidx(i)); %from the same group
    Adjacency_clusters(i,edges')=1;
end

G= digraph(f,t,edge_weights);
figure(1)
subplot(121)
h = plot(G);
title('Unnormalized Graph Laplacian and Non-weighted adjacency matrix ')
for i = 1:ncluster
    nodes = find(Qidx == i);
    [~,sl,al]= greedy_placement(QV,3,[],nodes',G,f,t,1);
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)
run_undersampled_optimization(G,casestr,sensor_locations,ammeter_locations);


%% Using the Y matrix Using the Real and imaginary Part of the Y matrix
close all;
G= digraph(f,t,edge_weights);
number_of_eigenvectors = 23;
U_k = U_Y_new (:,1:number_of_eigenvectors);
ncluster=8;
Q = full((Y_bus));
[QV, ~]=eigs(Q, ncluster,'smallestabs');
[Qidx, ~] = kmeans([real(QV) imag(QV)],ncluster);
%plot the adjacency:
ammeter_locations = [];
sensor_locations = [];
Adjacency_clusters=zeros(length(Qidx));
for i=1:length(Qidx)
    edges = find(Qidx==Qidx(i)); %from the same group
    Adjacency_clusters(i,edges')=1;
end

% [f,t] = findedge(G);
figure(2)
subplot(121)
h = plot(G);
title('Clustering by separating the real and imaginary of the Y matrix as Laplacian Matrix')
for i = 1:ncluster
    nodes = find(Qidx == i);
    [~,sl,al]= greedy_placement(U_k,3,[],nodes',G,f,t,1);
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
    
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


run_undersampled_optimization(G,casestr,sensor_locations,ammeter_locations);


%% Normalized laplacian and non-weighted adjacency matrix
G= graph(f,t,edge_weights);
degree_matrix = diag(degree(G,1:nb));
adjacency_matrix = adjacency(G);
laplacian_matrix_sym = speye(nb,nb) - (degree_matrix)^(-0.5)*adjacency_matrix*(degree_matrix)^(-0.5);
Q_sym = full(laplacian_matrix_sym);
[Q_L_sym, ~]=eigs(Q_sym, ncluster,'smallestabs');
Q_L_sym = Q_L_sym(:,1:ncluster);
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
figure(3)
subplot(121)
h = plot(G);
title('Clustering based on the Normalized Laplacian and non-weighted adjacency matrix')
for i = 1:ncluster
    nodes = find(Q_L_sym_idx == i);
    [~,sl,al]= greedy_placement(Q_sym,3,[],nodes',G,f,t,2);
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


run_undersampled_optimization(G,casestr,sensor_locations,ammeter_locations);



%% Laplacian using the covariance matrix

close all; 
PMU_cov = load(strcat(casestr,'_covariance','.mat')).covariance;
[PMUV, PMUD]  = eigs(PMU_cov, ncluster, 'largestabs');
% cluster as two different dimensions...
[PMUidx, PMUC,sumP] = kmeans([real(PMUV), imag(PMUV)], ncluster);
ammeter_locations = [];
sensor_locations = [];
Adjacency_clusters_cov=zeros(length(PMUidx));
for i=1:length(PMUidx)
  edges = find(PMUidx==PMUidx(i)); %from the same group
  Adjacency_clusters_cov(i,edges')=1;
end


figure(4)
subplot(121)
h = plot(G);
title('Clustering using the covariance matrix as the Laplacian matrix')
for i = 1:ncluster
    nodes = find(PMUidx == i);
    [~,sl,al]= greedy_placement(PMUV,3,[],nodes',G,f,t,1);
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


run_undersampled_optimization(G,casestr,sensor_locations,ammeter_locations);

%% Sensor Placement
number_of_eigenvectors = 23;
U_k = U_Y_new (:,1:number_of_eigenvectors);
% [idx,C] = kmeans(U_k,number_of_clusters);
[~,sensor_locations] = sensor_placement_matpower(mpc_case);
ammeter = load('ammeter.mat');

sensor_to_exclude=[76 56 43 22 85 77 22 47 11 13 80 69 40 45 38 67 37 23 66 84 82 78 15 11 54 55 51 18 20 71 75 16 39 34 36 42 74 79 72 59 11 62];
for s_e = 1:length(sensor_to_exclude)
    sensor_locations = sensor_locations(sensor_locations~=sensor_to_exclude(s_e));
end
sensor_locations = [sensor_locations [32 12 34 5 68 27 54 9 48 64 16 8 25 58 22 19 22 3]];



sensor_locations = sensor_locations';
sensor_locations = unique(sensor_locations);
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

c = Y_f*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb-1);
c = c + c_noise + 1i*c_noise;
ammeter_locations  = load('ammeter.mat').ammeter_location;
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
cm = abs(c(ammeter_locations));
M_a =zeros(length(ammeter_locations),nb-1);
for r = 1:length(ammeter_locations)
    M_a (r,ammeter_locations(r)) = 1; % Put 1 in sensor locations
end

% Select measurements according to sensor locations
sm= s(sensor_locations);
vm = abs(vcomplex(sensor_locations));
w_3 = 0;

figure(5)
h = plot(G);
title(strcat('For the best case: Number of Sensors:', string(number_of_sensors)))
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


%% Normalized laplacian and Weighted adjacency matrix
close all; 
ncluster=8;
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
h = plot(G);
title('Normalized Laplacian and Weighted adjacency matrix')
Yuk = Y_bus * U_Y_new;
% Yik=(F'*Y_f-T'*Y_t)./degree_vector'*U_Y_new;
Yik=(T'*Y_t+F'*Y_f)*U_Y_new;
n_eigs = 24;
for i = 1:ncluster
    nodes = find(Q_L_sym_idx == i);
    [~,sl,al]= greedy_placement(U_Y_new(:,1:n_eigs),2,[],nodes',G,f,t,3,Yuk(:,1:n_eigs),Yik(:,1:n_eigs));
    sensor_locations = [sensor_locations sl'];
    ammeter_locations = [ammeter_locations al];
    highlight(h,nodes,'NodeColor',C{i},'MarkerSize',6);
end
subplot(122)
h = plot(G);
highlight(h,sensor_locations,'NodeColor','k','MarkerSize',6);
highlight(h,f(ammeter_locations),t(ammeter_locations),'EdgeColor','red','LineWidth',1.5)


% 
% [~,sl,al]= greedy_placement(U_Y_new(:,1:n_eigs),64,[],1:nb,G,f,t,4,Yuk(:,1:n_eigs),Yik(:,1:n_eigs));
% G= graph(f,t,edge_weights);
% run_undersampled_optimization_buscurrent(G,casestr,sl);
G= graph(f,t,edge_weights);
run_undersampled_optimization_buscurrent(G,casestr,sensor_locations);


