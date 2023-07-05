% clc;
clear all;
% close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
ADD_NOISE = 1;
casestr = 'case85';
% casestr = 'case_ACTIVSg2000';
mpc_case = loadcase(casestr); % load case to get back to the original network

nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= digraph(f,t,edge_weights);

ncluster = 8;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.

G= graph(f,t,edge_weights);
F= full(sparse(1:G.numedges,f,1,G.numedges,G.numnodes));
T= full(sparse(1:G.numedges,t,1,G.numedges,G.numnodes));
degree_matrix = diag(degree(G,1:nb));
degree_vector = diag(degree_matrix);
degree_vector(degree_vector>1) = 2 ;
adjacency_matrix = adjacency(G);
adjacency_matrix = full(adjacency_matrix);
% degree_matrix = diag(sum(adjacency_matrix,2));
% laplacian_matrix = degree_matrix-full(adjacency_matrix);
laplacian_matrix_sym = speye(nb,nb) - (degree_matrix)^(-0.5)*adjacency_matrix*(degree_matrix)^(-0.5);
% laplacian_matrix_sym = speye(nb,nb) - (degree_matrix)^(-1)*adjacency_matrix;

% [Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
% mpcint = ext2int(mpc_case);
% Y_bus = makeYbus(mpcint);
% laplacian_matrix_sym = Y_bus;
Q_sym = full(laplacian_matrix_sym);

[U_Q_Y, Lambda_Q_Y] = eig(Q_sym);%./mpcint.baseMVA);
[val_Q_Y, ind_Q_Y]=sort(diag(Lambda_Q_Y),'ascend');
plot(abs(val_Q_Y))

rng('default');  % For reproducibility
eva = evalclusters(adjacency_matrix,'kmeans','silhouette','KList',[1:30])