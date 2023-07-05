clc;
clear all
close all 
define_constants;
casestr = 'case85';
mpc_case = loadcase(casestr); % load case to get back to the original network
mpc_case_results = runpf(mpc_case);
nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);
degree_vector = degree(G,1:G.numnodes);
degree_vector(degree_vector>1) = 2 ;
F= full(sparse(1:G.numedges,f,1,G.numedges,G.numnodes));
T= full(sparse(1:G.numedges,t,1,G.numedges,G.numnodes));
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180);
bus_current = (F'*Y_f-T'*Y_t)./degree_vector'*vcomplex;
% bus_current = bus_current./degree_vector';
% bus_injection = vcomplex.*conj(bus_current);
s = diag(vcomplex) * conj(Y_bus*vcomplex);
from_current = Y_f*vcomplex;
to_current = Y_t*vcomplex;

bus_to_current = T'*Y_t*vcomplex;
bus_current_2 = 





