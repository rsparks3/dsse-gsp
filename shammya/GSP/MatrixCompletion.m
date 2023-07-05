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
[C,ia,ib]=unique(f,'rows','stable');
repeat_entry = find(hist(ib,unique(ib))>1);

[C,ia,ib]=unique(t,'rows','stable');
repeat_entry_t = find(hist(ib,unique(ib))>1);

[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_bus = full(Y_bus);
[mpc_case_results,success] = runpf(mpc_case,opt);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180);
c_f = Y_f*vcomplex;
n_m = 12;
X = zeros(nbr,n_m);
X (:,1) = real(vcomplex(f));
X (:,2) = imag(vcomplex(f));

pd = mpc_case.bus(:,PD)/mpc_case.baseMVA;
qd = mpc_case.bus(:,QD)/mpc_case.baseMVA;

X (:,3) = real(c_f(f));
X (:,4) = imag(c_f(f));
X (:,5) = mpc_case_results.branch(:,PF);
X (:,6) = mpc_case_results.branch(:,QF);
X (:,7) = pd(f);
X (:,8) = qd(f);
X (:,9) = abs(vcomplex(f));
X (:,10) = real(vcomplex(t));
X (:,11) = imag(vcomplex(t));
X (:,12) = abs(vcomplex(t));

max_meas_avail = 12;
measurements = zeros(nbr,max_meas_avail);
for i = 1:nbr
    locs = randperm(n_m);
    measurements(i,:) = locs(1:max_meas_avail);
end
missing_row = 03;

X_meas = zeros(nbr-missing_row,max_meas_avail);
compare_matrix = zeros(nbr,n_m);
for i = 1:nbr-missing_row
    for j = 1:max_meas_avail
        X_meas(i,j) = X(i,measurements(i,j));
        compare_matrix(i,measurements(i,j)) = X(i,measurements(i,j));
    end
end

%%
cvx_clear
cvx_solver SDPT3
cvx_begin
    variable M(nbr,n_m)
    variable Y_meas(size(X_meas))
    variable epsilon_r(nbr,1) nonnegative
    variable epsilon_i(nbr,1) nonnegative
    variable gamma_f(nbr,1) nonnegative
    variable tau_p(nb,1) nonnegative
    variable tau_q(nb,1) nonnegative
    %     minimize (sum(epsilon) + sum(tau) + norm_nuc (M) )
    minimize (10*norm_nuc(M) + 2*sum(tau_p) + 2*sum(tau_q) + 2*sum(epsilon_r) + 2*sum(epsilon_i) + 2*sum(gamma_f))
    subject to
        
    for i = 1:nbr-missing_row
        for j = 1:max_meas_avail
            Y_meas(i,j) == M(i,measurements(i,j));
        end
    end
    norm(X_meas-Y_meas,'fro') <= 0.05;
    for i = 1:nb
        tau_p(i,1)  <= M(find(t==i),5) - sum(M(find(f==i),5)) - pd(i) <= tau_p(i,1);
        tau_q(i,1)  <= M(find(t==i),6) - sum(M(find(f==i),6)) - qd(i) <= tau_q(i,1);
    end
    for i = 1:nbr
        -epsilon_r(i,1) <= real(((M(i,1)+1i*M(i,2)) - (M(i,10)+1i*M(i,11))).*Y_bus(f(i),t(i)) - ((M(i,3)+1i*M(i,4)))) <= epsilon_r(i,1);
        -epsilon_i(i,1) <= imag(((M(i,1)+1i*M(i,2)) - (M(i,10)+1i*M(i,11))).*Y_bus(f(i),t(i)) - ((M(i,3)+1i*M(i,4)))) <= epsilon_i(i,1);
        -gamma_f <= M(i,12) - M(i,9) + mpc_case.branch(i,BR_R)* M(i,5) + mpc_case.branch(i,BR_X)* M(i,6) <= gamma_f(i,1);
    end
    for j = 1:length(repeat_entry)
        loc= find(f==repeat_entry(j));
        for k = 1:length(loc)-1
            M(loc(k),1) == M(loc(k+1),1);
            M(loc(k),2) == M(loc(k+1),2);
            M(loc(k),7) == M(loc(k+1),7);
            M(loc(k),8) == M(loc(k+1),8);
            M(loc(k),9) == M(loc(k+1),9);
        end
    end
    
    for j = 1:length(repeat_entry_t)
        loc= find(t==repeat_entry(j));
        for k = 1:length(loc)-1
            M(loc(k),10) == M(loc(k+1),10);
            M(loc(k),11) == M(loc(k+1),11);
            M(loc(k),12) == M(loc(k+1),12);
        end
    end

cvx_end
norm(M-X)
figure(1)
subplot(211)
plot(X(:,12),'linewidth',1.5)
hold on 
plot(M(:,12),'linewidth',1.5)
legend('Original','After Completion')
subplot(212)
plot(angle(X(:,10)+1i*X(:,11)),'linewidth',1.5)
hold on 
plot(angle(M(:,10)+1i*M(:,11)),'linewidth',1.5)
legend('Original','After Completion')
rms(X(:,12)-M(:,12))
rms(angle(X(:,10)+1i*X(:,11)) - angle(M(:,10)+1i*M(:,11)))


