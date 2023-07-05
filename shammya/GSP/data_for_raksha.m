clc;
clear all;
close all;
rng(1)   
define_constants;
INCLUDE_LOAD = 0;
casestr = 'case85.m';
mpc_case = loadcase(casestr);
base_case = mpc_case;
opt = mpoption('verbose',0,'out.all',0);
mpc_case_results = runpf(mpc_case);
number_of_scenarios = 5000;

V = zeros(length(mpc_case.bus),number_of_scenarios);
Vangle = zeros(length(mpc_case.bus),number_of_scenarios);
PD_case = zeros(length(mpc_case.bus),number_of_scenarios);
QD_case = zeros(length(mpc_case.bus),number_of_scenarios);

V(:,1) = mpc_case_results.bus(:,VM);
Vangle(:,1) =  mpc_case_results.bus(:,VA);
PD_case(:,1) = base_case.bus(:,PD);
QD_case(:,1) = base_case.bus(:,QD);
avg_pd = mean(base_case.bus(:,PD));
avg_qd = mean(base_case.bus(:,QD));

i = 1;
while i<number_of_scenarios
    mult = ones(1,length(mpc_case.bus));
%     Randomly choose the number of loads to be changed
    number_location = randi([1, length(mpc_case.bus)],1,1);
%     Randomly choose the locations of those loads
    locations = randperm(length(mpc_case.bus),number_location);
%     Get Multipliers between -0.5 and 1.5 and place them in the multiplier
%     matrix
    mult(1,locations) = (1.5+1.5)*rand(length(locations),1)-1.5;
%     Change real and reactive power 
    mpc_case.bus(:,PD) = base_case.bus(:,PD).*mult' ; 
    mpc_case.bus(:,QD) = base_case.bus(:,QD).*mult' ;
    for l = 1:length(locations)
        if mpc_case.bus(locations(l),PD) == 0
            mpc_case.bus(locations(l),PD) = avg_pd*mult(locations(l));
        end
        if mpc_case.bus(locations(l),QD) == 0
            mpc_case.bus(locations(l),QD) = avg_qd*mult(locations(l));
        end
        
    end
    
    [mpc_case_results,success] = runpf(mpc_case,opt);
    if success == 1
        fprintf('Completed Powr Flow For Scenario: %d\n',i);
        i = i+1;
        V(:,i) = mpc_case_results.bus(:,VM);
        Vangle(:,i) = mpc_case_results.bus(:,VA);
        PD_case(:,i) = mpc_case.bus(:,PD);
        QD_case(:,i) = mpc_case.bus(:,QD);
    end
end
j=1i;
% Formulate the complex voltage phasor
Vcomplex = zeros(length(mpc_case.bus),number_of_scenarios);
for r = 1:length(mpc_case.bus)
    for c = 1:number_of_scenarios
        Vcomplex (r,c) = V(r,c) * exp(j*Vangle(r,c)/180);
    end
    
end

nb = size(mpc_case.bus,1); %this is the number of buses
nl = size(mpc_case.branch,1); % this is the number of branches
% Calculate the co variance matrix if required
covaraince = Vcomplex * Vcomplex' / number_of_scenarios;


%% Convert to an orthogonal decomposition
[U_V_l,E_V] = eig(covaraince);
[val_V, ind_V]=sort(diag(E_V),'descend');
U_V_arr = U_V_l(:,ind_V');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_V_arr(:,i))*U_V_arr(:,i));
end
    
U_V_new = U_V_arr*diag(r_vect);

fprintf('Max Difference in estimate for Voltage:%f\n',max(abs(U_V_new * diag(val_V) * transpose(U_V_new) - covaraince),[],'all'));



[Y_bus,~,~] = makeYbus(mpc_case);
Y_L = mpc_case.bus(:,[PD,QD])./mpc_case.baseMVA;
Y_L = Y_L(:,1)-1i*Y_L(:,2);
Y_DD = full(Y_bus); % diagonal matrix 

[U_Y, Lambda_Y] = eig(Y_DD,'vector');%./mpcint.baseMVA);
[val_Y, ind_Y]=sort((Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(nb,1);
for i=1:nb
        r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end
    
U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));


%% Convert to Orthonormal decomposition
% The following process will not work as if number_of_eigenvectors = nb,
% the the GS outputs a U-K_perp that is an empty matrix
number_of_eigenvectors = 24;
U_k_gft = U_Y_new (:,1:number_of_eigenvectors);

% The orthonormal formulation using GS
[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);
U_k_perp = U_k_perp(:,1:number_of_eigenvectors);

V_tilde = zeros(number_of_eigenvectors,number_of_scenarios);
V_epsilon = zeros(number_of_eigenvectors,number_of_scenarios);
V_comb = zeros(length(mpc_case.bus),number_of_scenarios); % To save the sum of the voltages for further comparison
for i = 1:number_of_scenarios
    V_tilde(:,i) = (U_k)'*Vcomplex(:,i);
    V_epsilon(:,i) = (U_k_perp)'*Vcomplex(:,i);
    V_comb (:,i) = U_k * V_tilde(:,i) + U_k_perp*V_epsilon(:,i);
end

%%
x_axis_data = abs(val_Y)/max(abs(val_Y));
x_axis_data = x_axis_data(1:number_of_eigenvectors);
y_axis_mean = mean(abs(V_tilde),2);
y_axis_var = std(abs(V_tilde),0,2);
figure(1)
loglog(x_axis_data,y_axis_mean,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean+y_axis_var,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean-y_axis_var,'linewidth',1.5)
legend('mean','+variance','-variance')

figure(2)
y_axis_mean_ep = mean(abs(V_epsilon),2);
y_axis_var_ep = std(abs(V_epsilon),0,2);
loglog(x_axis_data,y_axis_mean_ep,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean_ep+y_axis_var_ep,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean_ep-y_axis_var_ep,'linewidth',1.5)
legend('mean','+variance','-variance')


% csvwrite('low_pass_85.csv',[x_axis_data y_axis_mean y_axis_mean+y_axis_var y_axis_mean-y_axis_var y_axis_mean_ep y_axis_mean_ep+y_axis_var_ep y_axis_mean_ep-y_axis_var_ep])


%% Three Phase Case
eigenvectors_reduction = 60;
BusSize = 34;
ADD_LOAD = 00;
[DSSCircObj, DSSText, gridpvPath] = DSSStartup;
DSSCircuit = DSSCircObj.ActiveCircuit;
DSSText.command = 'clear';
%%
[file_directory,slack_bus_name]= ChooseOpenDSSModel(BusSize);
DSSText.command = file_directory;
DSSText.Command = 'vsource.source.enabled=no';
DSSText.command = 'batchedit load..* enabled=no';
DSSText.Command = 'solve';

base_y = DSSCircuit.SystemY;
Ymatrix=reshape(base_y,[length(DSSCircuit.AllNodeNames)*2,length(DSSCircuit.AllNodeNames)]); % "2" is because the real/imaginary parts are separate
Ymatrix=Ymatrix';
Y_dss = Ymatrix(:,1:2:end) + sqrt(-1)*Ymatrix(:,2:2:end);

ynode_order = DSSCircuit.YNodeOrder;
[sorted_y_node,sort_index] =  sort(ynode_order);

sorted_y_dss= Y_dss(sort_index,:);
sorted_y_dss = sorted_y_dss(:,sort_index);

bus_names = cell(size(sorted_y_node));
for j = 1:length(sorted_y_node)
    name = split(sorted_y_node(j),'.');
    bus_names{j} = name{1};
end

bus_name_loc=containers.Map();
unique_bus_names = unique(bus_names,'stable');

for i = 1:length(unique_bus_names)
    bus_name_loc(unique_bus_names{i})=find(strcmpi(bus_names,unique_bus_names{i}));
end



%%
% Recompile the model to get a valid power flow
DSSText.command = file_directory;
DSSText.Command = 'solve';
% Retrieve all the voltages from the OpenDSS result
AllBusVoltage = DSSCircuit.AllbusVolts;
vcomplex = AllBusVoltage(1:2:end) + 1i * AllBusVoltage(2:2:end) ; % unit: Volt
vcomplex = transpose(vcomplex);
% sort the voltages the same way we sorted the y bus
sorted_vcomplex=vcomplex(sort_index);
% Calculate the original bus injection
s_original = diag(vcomplex) * conj(Y_dss*vcomplex)/1000;
% calculate the sorted bus injection values
s_sorted = diag(sorted_vcomplex) * conj(sorted_y_dss*sorted_vcomplex)/1000;
% Run a sanity check to ensure that the power calculations are correct
disp(max(abs(s_sorted-s_original(sort_index))));
Vcomplex = zeros(length(vcomplex),number_of_scenarios);
Base_Loads = getLoadInfo(DSSCircObj);
%%
for i = 1:number_of_scenarios
    number_location = randi([1, length(Base_Loads)],1,1);
    locations = randperm(length(Base_Loads),number_location);
    mult = (1.5+0.5)*rand(length(locations),1)-0.5;
    for l = 1:length(locations)
        DSSCircuit.Loads.Name = Base_Loads(l).name;
        DSSCircuit.Loads.kW = mult(l)*Base_Loads(l).kW;
        DSSCircuit.Loads.kvar = mult(l)*Base_Loads(l).kvar;
    end
    
    
    DSSText.Command = 'solve';
    AllBusVoltage = DSSCircuit.AllbusVolts;
    vcomplex = AllBusVoltage(1:2:end) + 1i * AllBusVoltage(2:2:end) ; % unit: Volt
    vcomplex = transpose(vcomplex);
% sort the voltages the same way we sorted the y bus
    sorted_vcomplex=vcomplex(sort_index);
    Vcomplex(:,i) = sorted_vcomplex;
    fprintf('Running iteration %d and minimum voltage:%f\n', i,abs(min(sorted_vcomplex)));
end

% Lets run the interpolation equation to see whether the sorted y_dss works or not
Y_DD = sorted_y_dss; % diagonal matrix

[U_Y, Lambda_Y] = eig(Y_DD);%./mpcint.baseMVA);
[val_Y, ind_Y]=sort(diag(Lambda_Y),'ascend');
U_Y_arr = U_Y(:,ind_Y');

% change into complex orthogonal matrices:
r_vect = zeros(length(Y_DD),1);
for i=1:length(Y_DD)
    r_vect(i) = 1./sqrt(transpose(U_Y_arr(:,i))*U_Y_arr(:,i));
end

U_Y_new = U_Y_arr*diag(r_vect);
fprintf('Max Difference in estimate for Y:%f\n',max(abs(U_Y_new * diag(val_Y) * transpose(U_Y_new) - Y_DD),[],'all'));

%% Phase Distribution
BusInfo = getBusInfo(DSSCircObj,unique_bus_names);
phase_distribution = zeros(length(unique_bus_names),3);
base_voltage = 1000*[BusInfo(:).kVBase]';

for i = 1:length(unique_bus_names)
    ph = BusInfo(i).nodes;
    mapping = bus_name_loc(upper(unique_bus_names{i}));
    phase_distribution(i,ph) = mapping;
end
%%
number_of_eigenvectors = length(Y_DD)-eigenvectors_reduction;
% vcomplex = full_voltage_magnitude;
U_k_gft = U_Y_new (:,1:number_of_eigenvectors) ;
[U_k,~,U_k_perp ] = gram_schemidt(U_k_gft);
U_k_perp = U_k_perp(:,:);
[~,c] = size(U_k_perp);
V_tilde = zeros(number_of_eigenvectors,number_of_scenarios);
V_epsilon = zeros(c,number_of_scenarios);
V_comb = zeros(length(vcomplex),number_of_scenarios); % To save the sum of the voltages for further comparison
for i = 1:number_of_scenarios
    V_tilde(:,i) = (U_k)'*Vcomplex(:,i);
    V_epsilon(:,i) = (U_k_perp)'*Vcomplex(:,i);
    V_comb (:,i) = U_k * V_tilde(:,i) + U_k_perp*V_epsilon(:,i);
end

figure()
plot(abs(Vcomplex(:,10)))
hold on 
plot(abs(V_comb(:,10)))

%%
x_axis_data = abs(val_Y)/max(abs(val_Y));
x_axis_data = x_axis_data(1:number_of_eigenvectors);
y_axis_mean = mean(abs(V_tilde),2);
y_axis_var = std(abs(V_tilde),0,2);
figure(3)
loglog(x_axis_data,y_axis_mean,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean+y_axis_var,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean-y_axis_var,'linewidth',1.5)
legend('mean','+variance','-variance')

figure(4)
y_axis_mean_ep = mean(abs(V_epsilon),2);
y_axis_var_ep = std(abs(V_epsilon),0,2);
loglog(x_axis_data,y_axis_mean_ep,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean_ep+y_axis_var_ep,'linewidth',1.5)
hold on
loglog(x_axis_data,y_axis_mean_ep-y_axis_var_ep,'linewidth',1.5)
legend('mean','+variance','-variance')