clc;
clear all;
close all;

define_constants;
ADD_NOISE = 0;
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_bus = full(Y_bus);
[mpc_case_results,success] = runpf(mpc_case,opt);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180);
c_f = Y_f*vcomplex;
n_m = 5;
X = zeros(nb,n_m);
X (:,1) = real(vcomplex);
X (:,2) = imag(vcomplex);
X (:,2) = abs(vcomplex);
pd = mpc_case.bus(:,PD)/mpc_case.baseMVA;
qd = mpc_case.bus(:,QD)/mpc_case.baseMVA;
X (:,4) = pd;
X (:,5) = qd;

max_meas_avail = 3;
measurements = zeros(nbr,max_meas_avail);
for i = 1:nb
    locs = randperm(n_m);
    measurements(i,:) = locs(1:max_meas_avail);
end

X_meas = zeros(nb,max_meas_avail);
compare_matrix = zeros(nbr,n_m);
for i = 1:nb
    for j = 1:max_meas_avail
        X_meas(i,j) = X(i,measurements(i,j));
        compare_matrix(i,measurements(i,j)) = X(i,measurements(i,j));
    end
end

%%
cvx_clear 

cvx_begin
    variable M(nb,n_m) 
    variable Y_meas(size(X_meas))
    variable tau_p(nb,1) nonnegative
    variable tau_q(nb,1) nonnegative
%     minimize (sum(epsilon) + sum(tau) + norm_nuc (M) )
    minimize (norm_nuc(M) + sum(tau_p) + sum(tau_p))
    subject to 
        norm(X_meas-Y_meas,'fro') <= 1;
        for i = 1:nb
            for j = 1:max_meas_avail
                Y_meas(i,j) == M(i,measurements(i,j));
            end
        end
        
%         for i = 1:nbr
%             -epsilon(i,1) <= abs(((M(f(i),1)+1i*M(f(i),2)) - (M(t(i),1)+1i*M(t(i),2))).*Y_bus(f(i),t(i)) - ((M(f(i),3)+1i*M(f(i),4)))) <= epsilon(i,1);
%            
%         end
% 
    for i = 1:nb
       tau_p(i,1)  <= M(find(t==i),5) - sum(M(find(f==i),6)) - M(i,4) <= tau_p(i,1);
       tau_q(i,1)  <= M(find(t==i),6) - sum(M(find(f==i),8)) - M(i,5) <= tau_q(i,1);
    end
cvx_end
 
