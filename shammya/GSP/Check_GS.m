clc;
clear all;
close all;
define_constants;
casestr = 'case22';
ADD_NOISE = 0;
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
mpc_case_results = runpf(mpc_case,opt); % run the power flow



%%
s_noise = ADD_NOISE*random_distribution(nb); % A customized function to generate random values
vcomplex = mpc_case_results.bus(:,VM).*exp(1j*mpc_case_results.bus(:,VA)*pi/180); % Preparing the complex voltages
s = diag(vcomplex) * conj(Y_bus*vcomplex)+s_noise+1i*s_noise; % Calculating power S = V* conjugate(Y*V) & add noise
v_noise = ADD_NOISE*random_distribution(nb);
v = vcomplex+v_noise+1i*v_noise;


number_of_eigenvectors = nb;
U_gft = U_Y_new (:,1:number_of_eigenvectors);

[U,S,~] = svd(U_gft);
[U_perp,~,~] = svd(eye(nb,nb)-U*U');

number_of_eigenvectors = nb-10;
U_k= U(:,1:number_of_eigenvectors);
U_k_perp= U_perp(:,1:number_of_eigenvectors);
v_tilde_1 = (U_k)'*v;
v_tilde_2 = (U_k_perp)'*v;

v_comb = (U_k * v_tilde_1 + U_k_perp*v_tilde_2);

figure(2)
subplot(211)
plot(abs(v_comb))
hold on
plot(abs(vcomplex))
legend('Reconstructed','Original')
subplot(212)
plot(angle(v_comb))
hold on
plot(angle(vcomplex))
legend('Reconstructed','Original')


