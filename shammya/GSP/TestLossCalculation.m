clc;
clear all;
close all;
rng(1)
define_constants;
bus_number = 100;
number_of_eigenvectors = 50;
gsp_approach = 1;

opt = mpoption('verbose',0,'out.all',0);
[n,e] = single_feeder_gen(bus_number,bus_number/20,0);

mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

generation_bus= randperm(bus_number,int16(bus_number/10));
% generation_bus (bus_number/10) = 100;

if sum(ismember(generation_bus,1))
else
   generation_bus(1) = 1 ;
end
generation_bus = sort(generation_bus,'ascend');

base_case = mpc_case;


nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
P_demand = mpc_case.bus(:,PD) / mpc_case.baseMVA;
Q_demand = mpc_case.bus(:,QD) / mpc_case.baseMVA;
mpc_case.bus(:,BS) = 0;

Q_shunt = mpc_case.bus(:,BS) / mpc_case.baseMVA;

G= graph(f,t);
mpc_case.branch(1,TAP) = 1;
% mpc_case.branch(:,BR_B) = 0;
[Y_bus,Y_f,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case,opt);
pgmax = sum(mpc_case.bus(:,PD))/mpc_case.baseMVA/length(generation_bus);
power_flow_limit = 1.5*mpc_case_results.branch(:,PF)/mpc_case.baseMVA;

vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise

sum((s_complex))
[loss,fchg,tchg ] = get_losses(mpc_case_results);
sum(loss)-1i*sum(fchg)-1i*sum(tchg)
% sum(mpc_case_results.bus(:,QD))
% mpc_case_results.gen(1,QG)

% mpc_case = base_case;
% mpc_case.bus(:,BS) = 0;
% mpc_case.branch(1,TAP) = 1;
% mpc_case.branch(:,BR_B) = 0;
% [Y_bus,Y_f,~] = makeYbus(mpc_case);
% [mpc_case_results,success] = runpf(mpc_case,opt);
% vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
% 
% s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise
% sum((s_complex))
% sum(get_losses(mpc_case_results))



