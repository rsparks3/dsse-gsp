clc;
clear all;
close all;
define_constants;
casestr = 'case85';
ADD_NOISE = 01;
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
base_case = mpc_case;
PD_case = base_case.bus(:,PD);
QD_case = base_case.bus(:,QD);

nb = length(mpc_case.bus);
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);

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



NUMBER_OF_EIGEN_SENSOR = 16;
U_k_gft_s = U_Y_new(:,1:NUMBER_OF_EIGEN_SENSOR);
[U_k_s,~, ] = gram_schemidt(U_k_gft_s);    
Yuk_s = Y_DD * U_k_s;

%%
mpc_case_results = runpf(mpc_case,opt); % run the power flow
s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;

c = Y_DD*vcomplex;
c_noise = ADD_NOISE*random_distribution(nb);
c = c + c_noise + 1i*c_noise;
vm_error = cell(3,1);
va_error = cell(3,1);
id = eye(nb,nb);
HV_base = [];

for i = 1:nb
    HV_base=[HV_base;transpose(kron(id(:,i),id(:,i)))];
end

HS_base = [];
for i = 1:nb
    HS_base=[HS_base;kron(transpose(id(:,i)),Y_DD(i,:))];
end

HC_base = [];
for i = 1:nb
    HC_base=[HC_base;kron(conj(Y_DD(i,:)),Y_DD(i,:))];
end
G_dir= digraph(f,t,edge_weights);
for bmult=0.25:0.25:0.75
    fprintf('\n Running Multiplier:%f\n',bmult);
    NUMBER_OF_SENSORS = ceil(nb*bmult);
    
    [~,sensor_locations,~]= greedy_placement(U_k_s,NUMBER_OF_SENSORS,[],1:nb,G_dir,f,t,1,Yuk_s,[]);
    cm = abs(c(sensor_locations));
    sm= conj(s(sensor_locations));
    vm = abs(v(sensor_locations));

    HV= HV_base(sensor_locations',:);
    HS= HS_base(sensor_locations',:);
    HC= HC_base(sensor_locations',:);
    
    xm=[vm.*vm;(sm);cm.*cm];
    H_M = [HV;HS;HC];
    avg_angle = containers.Map();avg_voltage=containers.Map();
    for number_of_eigenvectors = 10:5:nb
        
        U_k_gft = U_Y_new(:,1:number_of_eigenvectors);
        [U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);
        U_k_kron = kron(conj(U_k),U_k);
        H_tilde = H_M*U_k_kron;
        

        cvx_clear
            cvx_solver Mosek_3
            cvx_begin 
            cvx_precision high
            variable tildeW1(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
            minimize norm((xm-H_tilde*vec(tildeW1)))
            subject to
                H_tilde(1,:)*vec(tildeW1) == 1;
            cvx_end
        if strcmpi('Solved',cvx_status) || strcmpi('Inaccurate/Solved',cvx_status)
            [voltage_mag_new,voltage_angle_new] = calculate_voltage_phasor(tildeW1,U_k,G);
            angle_error = mean(abs(angle(vcomplex)-voltage_angle_new));
            voltage_error = mean(abs(abs(vcomplex) - voltage_mag_new));
            avg_angle(string(number_of_eigenvectors))  = angle_error ;
            avg_voltage(string(number_of_eigenvectors)) = voltage_error;

            fprintf('optval:%f RMS Angle:%f, RMS Voltage: %f  for Eigenvector:%d\n',cvx_optval,angle_error,voltage_error,number_of_eigenvectors);
        end
        pause(1)
    end
    
    vm_error{int8(bmult/0.25)} = avg_voltage;
    va_error{int8(bmult/0.25)} = avg_angle;
end
