clc;
clear all;
close all;
% rng(1)
define_constants;
bus_number = 400;
number_of_eigenvectors = 30;
slack_multi = 1;

opt = mpoption('verbose',0,'out.all',0);

[n,e] = single_feeder_gen(bus_number,bus_number/20,0);
mpc_case = matpower_fmt(n,e);
mpc_case = parallel_branch_join(mpc_case);

mpc_case = loadcase('case85.m');
bus_number = size(mpc_case.bus,1);

mpc_case.branch(:,BR_B) = 0;
mpc_case.bus(:,BS) = 0;
mpc_case.branch(1,TAP) = 1;
generation_bus= randperm(bus_number,int16(bus_number/10));
% generation_bus (bus_number/10) = 100;

if sum(ismember(generation_bus,1))
else
   generation_bus(1) = 1 ;
end
generation_bus = sort(generation_bus,'ascend');
% generation_bus=[1];




nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
P_demand = mpc_case.bus(:,PD) / mpc_case.baseMVA;
Q_demand = mpc_case.bus(:,QD) / mpc_case.baseMVA;
slack_bus = mpc_case.bus(:,BUS_TYPE)==3;

G= graph(f,t);
[Y_bus,~,~] = makeYbus(mpc_case);
[mpc_case_results,success] = runpf(mpc_case,opt);
pgmax = sum(mpc_case.bus(:,PD))/mpc_case.baseMVA/length(generation_bus);
power_flow_limit = 1.25*mpc_case_results.branch(:,PF)/mpc_case.baseMVA;

opf_mpc_case = mpc_case;
gen_row_base = opf_mpc_case.gen;
gencost_base = [2	0	0	3	0.11	0	0];
gencost = gencost_base;
gen_row = gen_row_base;
for r = 2:length(generation_bus)
    gen_row = [gen_row;gen_row_base];
    gencost = [gencost;gencost_base];
end


gen_row(:,1) = generation_bus';
gen_row (:,PMAX) = 10*pgmax*mpc_case.baseMVA;
gen_row (:,PMIN) = 0;
gen_row (1,PMAX) = 99;

opf_mpc_case.gen = gen_row;
b = 0.011; a = 0.85;
cost = (b-a).*rand(length(generation_bus),1) + a;
gencost(:,5) = cost;


opf_mpc_case.gencost = gencost;

%% Check All the Values
Y_DD = full(Y_bus);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
idx = sub2ind(size(Y_DD),mpc_case.branch(:,1),mpc_case.branch(:,2));
line_y = conj(-Y_DD(idx));
vv = vcomplex(mpc_case.branch(:,1))*(vcomplex(mpc_case.branch(:,1))-vcomplex(mpc_case.branch(:,2)))';
power_flow = mpc_case_results.baseMVA* diag(vv*diag((line_y)));

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise



%%
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
U_k_gft = U_Y_new (:,1:number_of_eigenvectors);


mult = zeros(nb,1);
mult(generation_bus) = 1;

% Swapping so that utility has the highest cost
Cs = cost;
[temp,loc] = max(Cs);
Cs(loc)= Cs(1);
Cs(1)= temp;

[U_k,~,~] = gram_schemidt(U_k_gft);

%% Covex Relaxation Code With GSP
A=U_k;
C=U_k';

H_f = [];
H_t = [];
H_tt = [];
for i = 1:length(f)
    H_f = [H_f;kron(transpose(C(:,f(i))),A(f(i),:))];
    H_t = [H_t;kron(transpose(C(:,t(i))),A(f(i),:))];
    H_tt = [H_tt;kron(transpose(C(:,t(i))),A(t(i),:))];
end
H_ft = H_f-H_t;

pgmax_vec = zeros(nb,1);
pgmax_vec(generation_bus) = 10*pgmax/100*randi([100 150],length(generation_bus),1);
id = speye(nb,nb);
[r_y,c_y] = size(Y_DD);
[r_i,c_i] = size (id);
H2 = zeros(nb,c_y*c_i);
for i = 1:nb
    H2(i,:)=kron(transpose(id(:,i)),Y_DD(i,:));
end
H2 = sparse(H2);

H1 = zeros(nb,c_i*c_i);

for i = 1:nb
    H1(i,:)=transpose(kron(id(:,i),id(:,i)));
end
H1 = sparse(H1);

U_k_kron = kron(conj(U_k),U_k);
H1_tilde = H1*U_k_kron;
H2_tilde = H2*U_k_kron;

disp('Matrix Multiplication Completed.')
%%

gsp_approach = 1;
if gsp_approach == 1
    tstart_gsp = tic;
    cvx_clear
    cvx_solver Mosek_3
    cvx_precision high
    cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable pg(nb) nonnegative
    variable qg(nb) 
    variable delta_p(nb)  
    variable delta_q(nb)
    minimize (sum(Cs.*pg(generation_bus).*pg(generation_bus))+1*norm(delta_p)+1*norm(delta_q) )
%     minimize (sum(Cs.*pg(generation_bus).*pg(generation_bus)))
    subject to
        0.85^2 <= H1_tilde*vec(tildeW) <= 1.1^2;
        pg(2:end) <= slack_multi*pgmax_vec(2:end);
        -0.75*pgmax_vec(2:end)*slack_multi <= qg(2:end) <= 0.75*pgmax_vec(2:end)*slack_multi;
        sum(pg-P_demand) - 1i*(sum(qg-Q_demand)) == sum(sum(H2_tilde*vec(tildeW)));
%         sum(delta_p)==0;
%         sum(delta_q)== 0;
%         sum(pg) >= sum(P_demand);
%         sum(qg) >= sum(Q_demand);
        (pg - P_demand) - 1i*(qg - Q_demand) - H2_tilde*vec(tildeW) == delta_p + 1i*delta_q;
        U_k(1,:)*(tildeW*U_k(1,:)') == 1.05^2;
%         -power_flow_limit <= real((H_ft*vec(tildeW)) .*line_y) <= power_flow_limit;
%         for br = 1:length(mpc_case.branch)
%             -power_flow_limit(br) <= real((U_k(f(br),:)*(tildeW*U_k(f(br),:)')-U_k(f(br),:)*(tildeW*U_k(t(br),:)')  )*(line_y(br))) <= power_flow_limit(br);
%         end
    cvx_end
    if strcmpi(cvx_status,'solved')
        telapsed_gsp = toc(tstart_gsp)
        [voltage_mag,voltage_angle] = calculate_voltage_phasor(tildeW,U_k,G);
        if slack_multi == 1
%         save(strcat('OPF_GSP_',string(bus_number),'.mat'),'bus_number','pg','qg','generation_bus','Cs','number_of_eigenvectors','mpc_case','voltage_mag','voltage_angle','telapsed_gsp','delta_p','delta_q')
        end
        
    end
end
%% re-solve the case
mpc_case_new_gsp=opf_mpc_case;
% mpc_case_new.gen(2:end,PG)= mpc_case_new.baseMVA*pg(generation_bus(2:end));
mpc_case_new_gsp.gen = mpc_case_new_gsp.gen(1,:);
mpc_case_new_gsp.gen(slack_bus,VG) = voltage_mag(1);
mpc_case_new_gsp.bus(generation_bus(2:end),PD) = mpc_case_new_gsp.bus(generation_bus(2:end),PD) - mpc_case_new_gsp.baseMVA*pg(generation_bus(2:end));
mpc_case_new_gsp.bus(generation_bus(2:end),QD) = mpc_case_new_gsp.bus(generation_bus(2:end),QD) - mpc_case_new_gsp.baseMVA*qg(generation_bus(2:end));
mpc_case_new_results_gsp = runpf(mpc_case_new_gsp,opt);
new_voltage = mpc_case_new_results_gsp.bus(:,VM);
[loss,fchg,tchg] = get_losses(mpc_case_new_results_gsp);
vcomplex = mpc_case_new_results_gsp.bus(:,VM).*exp(1j*mpc_case_new_results_gsp.bus(:,VA)*pi/180); % Preparing the complex voltages

s_complex = diag(vcomplex) * conj(Y_bus*vcomplex); % Calculating power S = V* conjugate(Y*V) & add noise

flow_gsp = real((H_ft*vec(tildeW)) .*line_y)*mpc_case_new_gsp.baseMVA;
flow_gsp_matpower = mpc_case_new_results_gsp.branch(:,PF);
fprintf('Voltage Magnitude Diff:%f\n',max(abs(voltage_mag-new_voltage)));
fprintf('Voltage Angle Diff:%f\n',max(abs(voltage_angle-angle(vcomplex))));
fprintf('Difference in Real Flow in MW:%f\n',max(abs( abs(flow_gsp)-abs(flow_gsp_matpower))))

figure(1)
subplot(311)
plot(new_voltage)
hold on 
plot(voltage_mag)
legend('matpower','voltage gsp')

subplot(312)
plot(angle(vcomplex))
hold on 
plot(voltage_angle)
legend('matpower','voltage gsp')

subplot(313)
plot(abs(flow_gsp_matpower))
hold on 
plot(abs(flow_gsp))
legend('matpower','voltage gsp')

