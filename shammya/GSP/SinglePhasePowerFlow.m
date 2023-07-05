clc;
clear all;
close all;
%% Define the global constants
SYN_FEEDER = false; % Whether the code is expecting a synthetic feeder or a general Matpower case, for general matpower  case, set to false, otherwise true 
opt = mpoption('verbose',0,'out.all',0);
define_constants; % Define the matpower constants
MATPOWER_CASE_FILE = 'case85'; % MATPOWER case file name, ignored if SYN_FEEDER = 1
BUS_NUMBER = 400; % Size of the synthetic feeder, ignored if SYN_FEEDER= 0
number_of_eigenvectors = 15; 
ORTHORNOMAL_DECOMPOSITION = 1; % 1 = Gram Schemidt, 2 = SVD
GSP = false; % if false, uses the basic convex relaxation equation and avoids GSP integration
cvx_solver SDPT3; % Choose appropriate solver for CVX
%% Load the appropriate case file
if SYN_FEEDER == true
    [n,e] = single_feeder_gen(BUS_NUMBER);
    mpc_case = matpower_fmt(n,e);
    mpc_case = parallel_branch_join(mpc_case);
    mpc_case.branch(:,BR_B) = 0; % Forcing the case BRANCH SUSCEPTENCE to 0
    mpc_case.bus(:,BS) = 0; % Shunt Reactance 0
    mpc_case.branch(1,TAP) = 1; % Force the TAP to be 1
else
    mpc_case = loadcase(MATPOWER_CASE_FILE);
end
% Retrieve required information for the case
nb = length(mpc_case.bus);
nbr = length(mpc_case.branch);
f = mpc_case.branch(:,F_BUS); % from nodes
t = mpc_case.branch(:,T_BUS); % to nodes
G= graph(f,t); % The graph will be required for the voltage calculation part
P_demand = mpc_case.bus(:,PD) / mpc_case.baseMVA; % real power demand
Q_demand = mpc_case.bus(:,QD) / mpc_case.baseMVA; % reactive power demand
slack_bus = find(mpc_case.bus(:,BUS_TYPE)==3); % locate the slack bus
non_slack_bus = find(mpc_case.bus(:,BUS_TYPE)~=3);
[mpc_case_results,success] = runpf(mpc_case,opt); % Run the power flow

if success ~= 1
    error('Power Flow Did not converge, check your network model.');
end
[Y_bus,Y_f,Y_t] = makeYbus(mpc_case);
Y_DD = full(Y_bus);
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s_complex = diag(vcomplex) * conj(Y_bus*vcomplex);
load_buses = find(mpc_case.bus(:,PD));

slack_bus_voltage = abs(vcomplex(slack_bus));
idx = sub2ind(size(Y_DD),mpc_case.branch(:,1),mpc_case.branch(:,2));
line_y = conj(-Y_DD(idx));
vv = vcomplex(mpc_case.branch(:,1))*(vcomplex(mpc_case.branch(:,1))-vcomplex(mpc_case.branch(:,2)))';
power_flow = mpc_case_results.baseMVA* diag(vv*diag((line_y))); % calculate the power flow for lines using voltage difference equation

%% Retrieve information for the Graph Shift Operator
[U_Y, Lambda_Y] = eig(Y_DD);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices
r_vect = zeros(nb,1);
for i=1:nb
    r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end

U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all')); % This difference should be very low
if number_of_eigenvectors > length(U_Y_new)
    error('Number of Eigenvectors > Size of Y matrix')
end
U_k_gft = U_Y_new (:,1:number_of_eigenvectors);
if ORTHORNOMAL_DECOMPOSITION  == 1
    [U_k,~,~] = gram_schemidt(U_k_gft);
elseif ORTHORNOMAL_DECOMPOSITION  == 2
    [U_k,~,~] = svd(U_k_gft);
    U_k = U_k (:,1:number_of_eigenvectors);
else
    error('Orthonormal Decomposition Mode not Understood, Use 1 for Gram Schemidt, 2 for SVD')
end

%% 
GSP = false;
if GSP == false
    cvx_clear
    cvx_precision low
    cvx_begin 
       variable W(nb,nb) complex semidefinite
       variable pg_convex(nb) nonnegative
       variable qg_convex(nb) 
       minimize pg_convex(slack_bus)
       subject to
           (pg_convex-P_demand) - 1i* (qg_convex-Q_demand) == diag(Y_DD*W); % The negative sign appears to account for the conjugate on the RHS
           pg_convex(non_slack_bus) ==0; % non slack buses have no generation
           qg_convex(non_slack_bus) ==0; % non slack buses have no generation
           W(1,1) == slack_bus_voltage^2; % Force the slack bus constraint
    cvx_end

    if strcmpi(cvx_status,'solved')
            [voltage_mag_gsp,voltage_angle_convex] = calculate_voltage_phasor(W,U_k,G);
            mae_angle_gsp  = mae(angle(vcomplex)-voltage_angle_convex) ;
            mae_voltage = mae(abs(vcomplex) - voltage_mag_gsp);
            fprintf('optval:%f MAE Angle:%f, MAE Voltage: %f\n',cvx_optval,mae_angle_gsp,mae_voltage);
            figure(1)
            subplot(211)
            plot(voltage_mag_gsp)
            hold on
            plot(abs(vcomplex))
            legend('Reconstructed','Original')
            subplot(212)
            plot(voltage_angle_convex)
            hold on
            plot(angle(vcomplex))
            legend('Reconstructed','Original')
    else
       error('CVX Could not Solve the convex optimization problem.'); 
    end
end

%%

if GSP == true
    id = speye(nb,nb);
    [r_y,c_y] = size(Y_DD);
    [r_i,c_i] = size (id);
    H1 = zeros(nb,c_i*c_i);
    s_demand = conj(P_demand + 1i*Q_demand);
    for i = 1:nb
        H1(i,:)=transpose(kron(id(:,i),id(:,i)));
    end
    H1 = sparse(H1);
    
    H2 = zeros(nb,c_y*c_i);
    for i = 1:nb
        H2(i,:)=kron(transpose(id(:,i)),Y_DD(i,:));
    end
    H2 = sparse(H2);

    U_k_kron = kron(conj(U_k),U_k);
    H1_tilde = H1*U_k_kron;
    H2_tilde = H2*U_k_kron;
    cvx_clear
    cvx_solver SDPT3
    cvx_begin 
       variable W_gsp(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
       variable pg_gsp(nb) nonnegative
       variable qg_gsp(nb) 
       variable delta_p(nb)  
       variable delta_q(nb) 
       minimize ( pg_gsp(slack_bus) + 100*norm(delta_p) + 100*norm(delta_q))
       
       subject to
            0.95^2 <= H1_tilde*vec(W_gsp) <= 1.05^2;
           (pg_gsp - P_demand) - 1i*(qg_gsp - Q_demand) - H2_tilde*vec(W_gsp) == delta_p + 1i*delta_q; % The negative sign appears to account for the conjugate on the RHS
           pg_gsp(non_slack_bus) ==0; 
           qg_gsp(non_slack_bus) ==0;
           sum(pg_gsp) >= sum(P_demand);
           sum(qg_gsp) >= sum(Q_demand);
           U_k(1,:)*(W_gsp*U_k(1,:)') == slack_bus_voltage^2;  

    cvx_end
    if strcmpi(cvx_status,'solved')
            [voltage_mag_gsp,voltage_angle_gsp] = calculate_voltage_phasor(W_gsp,U_k,G);
            mae_angle_gsp  = mae(angle(vcomplex)-voltage_angle_gsp) ;
            mae_voltage_gsp = mae(abs(vcomplex) - voltage_mag_gsp);
            fprintf('optval:%f MAE Angle:%f, MAE Voltage: %f\n',cvx_optval,mae_angle_gsp,mae_voltage_gsp);
            figure(2)
            subplot(211)
            plot(voltage_mag_gsp)
            hold on
            plot(abs(vcomplex))
            legend('Reconstructed','Original')
            subplot(212)
            plot(voltage_angle_gsp)
            hold on
            plot(angle(vcomplex))
            legend('Reconstructed','Original')
    else
       error('CVX Could not Solve the convex optimization problem.'); 
    end
end


%%
Uk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors^2);
YUk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors^2);
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
s_demand = (P_demand +1i*Q_demand);
cvx_clear
cvx_begin 
    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb-1) complex
    minimize norm(s_est-s_demand(2:end))
    for k = 2:nb
        s_est(k-1) == conj(YUk_UkH(k,:)*vec(transpose(tildeW))); 
    end
    U_k(1,:)*(tildeW*U_k(1,:)') == slack_bus_voltage^2; 
cvx_end
if strcmpi(cvx_status,'solved')
    [voltage_mag_gsp,voltage_angle_gsp] = calculate_voltage_phasor(tildeW,U_k,G);
    mae_angle_gsp  = mae(angle(vcomplex)-voltage_angle_gsp) ;
    mae_voltage_gsp = mae(abs(vcomplex) - voltage_mag_gsp);
    fprintf('optval:%f MAE Angle:%f, MAE Voltage: %f\n',cvx_optval,mae_angle_gsp,mae_voltage_gsp);
    figure(3)
    subplot(211)
    plot(voltage_mag_gsp)
    hold on
    plot(abs(vcomplex))
    legend('Reconstructed','Original')
    subplot(212)
    plot(voltage_angle_gsp)
    hold on
    plot(angle(vcomplex))
    legend('Reconstructed','Original')
else
    error('CVX Could not Solve the convex optimization problem.'); 
end

