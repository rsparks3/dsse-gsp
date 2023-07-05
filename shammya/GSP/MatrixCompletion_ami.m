clc;
clear all;
close all;

define_constants;
ADD_NOISE = 0;
casestr = 'case141';
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

Y11 = Y_bus(1,1);
Y1L = Y_bus(1,2:end);
YL1 = Y_bus(2:end,1);
YLL = Y_bus(2:end,2:end);
s = diag(vcomplex) * conj(Y_bus*vcomplex);
s1 = s(1);
v1 = vcomplex(1);
v_1 = vcomplex(2:end);
w = -v1*inv(YLL)*YL1;
A = [inv(YLL)*inv(diag(conj(w))) -1i*inv(YLL)*inv(diag(conj(w)))];
C = inv(diag(abs(w))) * real(diag(conj(w))*A);

% tau_r = real(v_1-A*[real(s(2:end)); imag(s(2:end))]-w);
% tau_i = imag(v_1-A*[real(s(2:end)); imag(s(2:end))]-w);
% gamma = abs(v_1) - ( C*[real(s(2:end)); imag(s(2:end))] + abs(w));
% alpha_r = real(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(v_1))));
% alpha_i = imag(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(v_1))));


c = Y_bus*vcomplex;
c = abs(c);
c_f = Y_f*vcomplex;
n_m = 5;
% ami_buses = find(mpc_case.bus(:,PD)>0);

ami_buses = randperm(nb,ceil(0.5*nb));
% ami_buses (1)=83;
nonami_buses = setdiff(1:nb,ami_buses);
X = zeros(nb,n_m);
X (:,1) = real(vcomplex);
X (:,2) = imag(vcomplex);

pd = mpc_case.bus(:,PD)/mpc_case.baseMVA;
qd = mpc_case.bus(:,QD)/mpc_case.baseMVA;


X (:,4) = real(s);
X (:,5) = imag(s);
X (:,3) = abs(vcomplex);
% X (:,9) = abs(vcomplex(f));
% X (:,10) = real(vcomplex(t));
% X (:,11) = imag(vcomplex(t));
% X (:,12) = abs(vcomplex(t));

max_meas_avail = 5;
% measurements = zeros(nb,max_meas_avail);
% for i = 1:nb
%     locs = randperm(n_m);
%     measurements(i,:) = locs(1:max_meas_avail);
% end

measurements = X;


% X_meas = zeros(length(ami_buses),max_meas_avail);

X_meas = measurements(ami_buses,:);

abs_w = abs(w);
%%
clear v_1
cvx_clear
cvx_solver SDPT3
cvx_begin
    variable M(size(X))
    variable Y_meas(size(X_meas))
    variable tau_r(nb-1,1) nonnegative
    variable tau_i(nb-1,1) nonnegative
    variable alpha_r nonnegative
    variable alpha_i nonnegative
    variable v_1(nb-1) complex
    variable gamma_f(nb-1) nonnegative
    %     minimize (sum(epsilon) + sum(tau) + norm_nuc (M) )
    minimize (norm_nuc(M) + 10*norm(tau_r) + 10*norm(tau_i) + 10*norm(gamma_f) + 10*norm(alpha_r) + 10*norm(alpha_i))
    subject to
        norm(X_meas-M(ami_buses,:),'fro') <= 0.0001;

        M(1,1) +1i*M(1,2) == v1; 

        v_1 == M(2:end,1) + 1i*M(2:end,2);
        -tau_r <= real(v_1-A*[M(2:end,4); M(2:end,5)]-w)  <= tau_r;
        -tau_i <= imag(v_1-A*[M(2:end,4); M(2:end,5)]-w)  <= tau_i;
        -gamma_f <= (M(2:end,3) - ( C*[M(2:end,4); M(2:end,5)] + abs_w)) <= gamma_f;
        -alpha_r<= real(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(v_1)))) <= alpha_r;
        -alpha_i<= imag(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(v_1)))) <= alpha_i;
        
%         -tau_r <= real(M(2:end,1) + 1i*M(2:end,2)-A*[M(2:end,4); M(2:end,5)]-w)  <= tau_r;
%         -tau_i <= imag(M(2:end,1) + 1i*M(2:end,2)-A*[M(2:end,4); M(2:end,5)]-w)  <= tau_i;
%         -gamma_f <= (M(2:end,3) - ( C*[M(2:end,4); M(2:end,5)] + abs_w)) <= gamma_f;
%         -alpha_r<= real(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(M(2:end,1) + 1i*M(2:end,2))))) <= alpha_r;
%         -alpha_i<= imag(s1-(v1*(conj(Y11)*conj(v1)+conj(Y1L)*conj(M(2:end,1) + 1i*M(2:end,2))))) <= alpha_i;
    

%     for i = 1:nb
%         tau_p(i,1)  <= M(find(t==i),5) - sum(M(find(f==i),5)) - pd(i) <= tau_p(i,1);
%         tau_q(i,1)  <= M(find(t==i),6) - sum(M(find(f==i),6)) - qd(i) <= tau_q(i,1);
%     end
%     for i = 1:nbr
%         -epsilon_r(i,1) <= real(((M(i,1)+1i*M(i,2)) - (M(i,10)+1i*M(i,11))).*Y_bus(f(i),t(i)) - ((M(i,3)+1i*M(i,4)))) <= epsilon_r(i,1);
%         -epsilon_i(i,1) <= imag(((M(i,1)+1i*M(i,2)) - (M(i,10)+1i*M(i,11))).*Y_bus(f(i),t(i)) - ((M(i,3)+1i*M(i,4)))) <= epsilon_i(i,1);
%         -gamma_f <= M(i,12) - M(i,9) + mpc_case.branch(i,BR_R)* M(i,5) + mpc_case.branch(i,BR_X)* M(i,6) <= gamma_f(i,1);
%     end
%     for j = 1:length(repeat_entry)
%         loc= find(f==repeat_entry(j));
%         for k = 1:length(loc)-1
%             M(loc(k),1) == M(loc(k+1),1);
%             M(loc(k),2) == M(loc(k+1),2);
%             M(loc(k),7) == M(loc(k+1),7);
%             M(loc(k),8) == M(loc(k+1),8);
%             M(loc(k),9) == M(loc(k+1),9);
%         end
%     end
%     
%     for j = 1:length(repeat_entry_t)
%         loc= find(t==repeat_entry(j));
%         for k = 1:length(loc)-1
%             M(loc(k),10) == M(loc(k+1),10);
%             M(loc(k),11) == M(loc(k+1),11);
%             M(loc(k),12) == M(loc(k+1),12);
%         end
%     end

cvx_end
norm(M-X)
figure(1)
subplot(211)
plot(abs(X(:,1)+1i*X(:,2)),'linewidth',1.5)
hold on 
plot(abs(M(:,1)+1i*M(:,2)),'linewidth',1.5)
legend('Original','After Completion')
subplot(212)
plot(angle(X(:,1)+1i*X(:,2)),'linewidth',1.5)
hold on 
plot(angle(M(:,1)+1i*M(:,2)),'linewidth',1.5)
legend('Original','After Completion')
% rms(X(:,12)-M(:,12))
% rms(angle(X(:,10)+1i*X(:,11)) - angle(M(:,10)+1i*M(:,11)))


