clc;
clear all;
close all;
eigenvectors_reduction = 40;
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
U_k = U_Y_new (:,1:number_of_eigenvectors) ;
v_new = U_k*transpose(U_k)*sorted_vcomplex;
figure(3)
% subplot(311)
% loc = find(abs(full_voltage_magnitude(loc))>0);
subplot(211)
plot(abs(v_new),'linewidth',1.5)
hold on
plot(abs(sorted_vcomplex),'linewidth',1.5)
legend('Interpolated','Original')
subplot(212)
plot(angle(v_new),'linewidth',1.5)
hold on
plot(angle(sorted_vcomplex),'linewidth',1.5)
legend('Interpolated','Original')




%% Sensor placement

% Yuk = sorted_y_dss * U_Y_new;
% Yuk = Yuk(:,1:number_of_eigenvectors);
% M_s = [];
% M= 52;
% M_tilde = M-length(M_s); % total number to be found
% [N,~]=size(U_k);
% i=1;
% j_tilde = 1:length(sorted_y_dss);
% j_tilde(M_s)=[];
% while i <= M_tilde
%     sigma_min_uk=zeros(length(j_tilde),1);
%     sigma_min_yuk=zeros(length(j_tilde),1);
%     for j=1: length(j_tilde)
%         A = U_k([M_s;j_tilde(j)],:);
%         [~,sigma_min_uk(j),~,flag] = svds(A,1,'smallest');
%         if(flag)
%             continue;
%         end
%         B = Yuk([M_s;j_tilde(j)],:);
%         [~,sigma_min_yuk(j),~,flag] = svds(B,1,'smallest');
%         if(flag)
%             continue;
%         end
%     end
% %         [~,min_ind]=max(abs(sigma_min_uk)+abs(sigma_min_yuk));
%     [~,min_ind]= max(min([sigma_min_uk sigma_min_yuk],[],2));
%     M_s = [M_s;j_tilde(min_ind)];
%     j_tilde(min_ind)=[];
%     i=i+1;
% end
% P_S = eye(N);
% P_S=P_S(:,M_s');
%
% bus_name_sensor = unique(bus_names(M_s));
%
% sensor_locations = zeros(length(bus_names),1);
%
% for i = 1:length(bus_name_sensor)
%     x = find(strcmp(bus_names,bus_name_sensor(i)));
%     sensor_locations(x,1)=1;
% end


%% Sensor placement Bus Wise



%% Construct the graph for node relabeling
from_bus_nodes = [];
to_bus_nodes = [];
lines = getLineInfo(DSSCircObj);
transformers = getTransformerInfo(DSSCircObj);
bus_names = cell(length(lines)*2 + length(transformers)*2,1);
counter = 1;
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    
    bus_names{counter} = name{1};
    counter = counter +1 ;
    name = split(lines(i).bus2,'.');
    bus_names{counter} =name{1};
    counter = counter +1 ;
end

for i = 1:length(transformers)
    name = split(transformers(i).bus1,'.');
    bus_names{counter} = name{1};
    counter = counter +1 ;
    name = split(transformers(i).bus2,'.');
    bus_names{counter} = name{1};
    counter = counter +1 ;
end


bus_info=getBusInfo(DSSCircObj,unique_bus_names);
inner_map = containers.Map(unique_bus_names,1:length(unique_bus_names));
% looping through the lines for finding all the from and nodes
for i = 1:length(lines)
    name= split(lines(i).bus1,'.');
    from_bus_loc = inner_map(upper(name{1}));
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(lines(i).bus2,'.');
    to_bus_loc = inner_map(upper(name{1}));
    to_bus_nodes = [to_bus_nodes to_bus_loc];
end
% Append the from and to nodes of transformers to the corresponding list
for i = 1:length(transformers)
    name= split(transformers(i).bus1,'.');
    from_bus_loc = inner_map(upper(name{1}));
    from_bus_nodes = [from_bus_nodes from_bus_loc];
    name= split(transformers(i).bus2,'.');
    to_bus_loc = inner_map(upper(name{1}));
    to_bus_nodes = [to_bus_nodes to_bus_loc];
end
%% Process the graph
remove_list=[];
for i = 1:length(from_bus_nodes)
    if ~ismember(i,remove_list)
        f=from_bus_nodes(i);
        t= to_bus_nodes(i);
        for j = 1:length(from_bus_nodes)
            if j ~=i
                if (from_bus_nodes(j)==f && to_bus_nodes(j)==t)
                    remove_list = [remove_list j];
                end
            end
        end
    end
end

from_bus_nodes(remove_list) = [];
to_bus_nodes(remove_list) = [];

[nmap, rmap, fnew, tnew] = NodeRelabling(from_bus_nodes',to_bus_nodes',inner_map(upper(slack_bus_name)));
G_old = graph(from_bus_nodes,to_bus_nodes);




%%
% node_ordering_opendss = DSSCircuit.YNodeOrder;
% node_ordering_opendss_names = cell(size(node_ordering_opendss));
G = digraph(fnew,tnew);









% for i = 1:length(node_ordering_opendss)
%     name= split(node_ordering_opendss(i),'.');
%     node_ordering_opendss_names{i} = name{1};
% end
% opendss_bus_names = DSSCircuit.AllBusNames;
testkeys = keys(inner_map);

edge_weights = [];
%% Clustering of nodes
for i = 1:length(fnew)
    f = fnew(i); t = tnew(i);
    fp=rmap(f); tp = rmap(t);
    testind = cellfun(@(x)isequal(x,fp),values(inner_map));
    fpp = upper(testkeys(testind));
    testind = cellfun(@(x)isequal(x,tp),values(inner_map));
    tpp = upper(testkeys(testind));
    from_bus_loc = bus_name_loc(fpp{1});
    to_bus_loc = bus_name_loc(tpp{1});
    line_y = sorted_y_dss(from_bus_loc,:);
    line_y = line_y(:,to_bus_loc');
    [~,c]= size(line_y);
    if c > 1
        edge_weights = [edge_weights;mean(abs(diag(line_y)))];
    else
        edge_weights = [edge_weights;mean(abs((line_y)))];
    end
end
G_weighted = graph(fnew,tnew);

F= full(sparse(1:G_weighted.numedges,f,1,G_weighted.numedges,G_weighted.numnodes));
T= full(sparse(1:G_weighted.numedges,t,1,G_weighted.numedges,G_weighted.numnodes));
degree_matrix = diag(degree(G_weighted));
degree_vector = diag(degree_matrix);
degree_vector(degree_vector>1) = 2 ;
adjacency_matrix = adjacency(G_weighted,'weighted');
% degree_matrix = diag(sum(adjacency_matrix,2));
% laplacian_matrix = degree_matrix-full(adjacency_matrix);

laplacian_matrix_sym = speye(G_weighted.numnodes,G_weighted.numnodes) - (degree_matrix)^(-0.5)*adjacency_matrix*(degree_matrix)^(-0.5);
Q_sym = full(laplacian_matrix_sym);
nb = G_weighted.numnodes;
ncluster = 6;
sensors_per_cluster = 4;
[Q_L_sym_whole, ~]=eigs(Q_sym, ncluster,'smallestabs');
Q_L_sym = Q_L_sym_whole;
% ammeter_locations = [];
% sensor_locations = [];
for i = 1:nb
    Q_L_sym(i,:) =  Q_L_sym(i,:) / sqrt(sum(Q_L_sym(i,:).*Q_L_sym(i,:)));
end
[Q_L_sym_idx, ~] = kmeans(Q_L_sym,ncluster);

Adjacency_clusters_L_sym=zeros(length(Q_L_sym_idx));
for i=1:length(Q_L_sym_idx)
    edges = find(Q_L_sym_idx==Q_L_sym_idx(i)); %from the same group
    Adjacency_clusters_L_sym(i,edges')=1;
end
G= digraph(fnew,tnew,edge_weights);
figure(7)
subplot(121)
h = plot(G);
title('Normalized Laplacian and Weighted adjacency matrix')
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6],'m',[.8 .2 .9],[.1 .4 .7]}; % Cell array of colors.
Yuk = sorted_y_dss * U_Y_new;
Yuk = Yuk(:,1:number_of_eigenvectors);
M_ss = [];
M_s_y_s = [];
for ii = 1:ncluster
    nodes = find(Q_L_sym_idx == ii);
    original_node_map = rmap(nodes);
    sensor_unique_bus_names = cell(length(original_node_map),1);
    for nn = 1:length(original_node_map)
        testind = cellfun(@(x)isequal(x,original_node_map(nn)),values(inner_map));
        sensor_unique_bus_names(nn,1) = upper(testkeys(testind));
    end
    M= sensors_per_cluster;
    M_s = [];
    M_s_y= [];
    M_tilde = M-length(M_s); % total number to be found
    [N,~]=size(U_k);
    i=1;
    while i <= M_tilde
        sigma_min_uk=zeros(length(sensor_unique_bus_names),1);
        sigma_min_yuk=zeros(length(sensor_unique_bus_names),1);
        sigma_min_uk(M_s) = nan; sigma_min_yuk(M_s) = nan;
        for j=1: length(sensor_unique_bus_names)
            if ~any(M_s == original_node_map(j))
                sensor_bus_map = bus_name_loc(sensor_unique_bus_names{j});
                A = U_k([M_s_y sensor_bus_map'],:);
                [~,sigma_min_uk(j),~,flag] = svds(A,1,'smallest');
                if(flag)
                    continue;
                end
                B = Yuk([M_s_y sensor_bus_map'],:);
                [~,sigma_min_yuk(j),~,flag] = svds(B,1,'smallest');
                if(flag)
                    continue;
                end
            end
        end
        [~,min_ind]= max(min([sigma_min_uk sigma_min_yuk],[],2));
        M_s = [M_s;original_node_map(min_ind)];
        
        i=i+1;
    end
    %     P_S = eye(N);
    %     P_S=P_S(:,M_s');
    M_s_y = [];
    for k = 1:length(M_s)
        sensor_placed_map = bus_name_loc(unique_bus_names{M_s(k)});
        M_s_y = [M_s_y  sensor_placed_map'];
    end
    M_ss = [M_ss; M_s];
    M_s_y_s = [M_s_y_s M_s_y];
    highlight(h,nodes,'NodeColor',C{ii},'MarkerSize',6);
end
bus_name_sensor = unique_bus_names(M_ss);
total_sensors = length(M_ss);
sensor_locations = M_s_y_s;
graph_highlight = zeros(1,length(bus_name_sensor));
for i = 1:length(graph_highlight)
    name = bus_name_sensor(i);
    graph_highlight(1,i) = nmap(inner_map(upper(name{1})));
end





%%
nb = length(sorted_y_dss);
% number_of_sensors = sum(sensor_locations);
% sensor_locations= find(sensor_locations==1);
number_of_sensors = length(sensor_locations);
% sensor_locations = 1:nb;
M =zeros(number_of_sensors,nb);
for r = 1:length(sensor_locations)
    M (r,sensor_locations(r)) = 1; % Put 1 in sensor locations
end
subplot(122)
h= plot(G);
highlight(h,graph_highlight,'Nodecolor','r');
%% Covex Relaxation Code
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

C = U_k'*Y_DD';
for k = 1:nb
   UkH_YkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end
% clear tildeW
%%
sm = s_sorted*1000;
sm = sm(sensor_locations);
vm = abs(sorted_vcomplex);
vm=vm(sensor_locations);
c = sorted_y_dss*sorted_vcomplex;
cm = abs(c);
cm = cm(sensor_locations);
cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin

variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
variable s_est(nb) complex
variable v_est(nb)
variable a_est(nb)
minimize (25*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

subject to

    tildeW == tildeW';
%     s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%     v_est == diag(U_k*tildeW*(U_k)');
    for k = 1:nb
        v_est(k) == Uk_UkH(k,:)*vec(transpose(tildeW));
        s_est(k) == conj( YUk_UkH(k,:)*vec(transpose(tildeW)) );
        a_est(k) == UkH_YkH(k,:)*vec(transpose(tildeW)); 
    end
%     0 <= sum(real(s_est)) <=2e+05;
%     0 <= sum(imag(s_est)) <=5e+05;
%     a_est == diag(Y_DD*U_k*tildeW*U_k'*Y_DD');
    U_k(bus_name_loc(slack_bus_name),:)*(tildeW*U_k(bus_name_loc(slack_bus_name),:)') == sorted_vcomplex(bus_name_loc(slack_bus_name))*sorted_vcomplex(bus_name_loc(slack_bus_name))';
cvx_end
W = U_k * tildeW *(U_k)';



%% Voltage Reconstruction from W
if strcmpi('solved',cvx_status)
    
    all_voltage = zeros(size(sorted_vcomplex));
    all_voltage(bus_name_loc(slack_bus_name),1) = sorted_vcomplex(bus_name_loc(slack_bus_name));
    for node = 2:length(rmap)
        to = node;
        from = predecessors(G,to);
        to = (find(nmap==to));
        testind = cellfun(@(x)isequal(x,to),values(inner_map));
        to_msg_key = upper(testkeys(testind));
        to= find(strcmpi(unique_bus_names,to_msg_key));
        
        from = (find(nmap==from));
        testind = cellfun(@(x)isequal(x,from),values(inner_map));
        from_msg_key = upper(testkeys(testind));
        
        from= find(strcmpi(unique_bus_names,from_msg_key));
        source_voltage = all_voltage(bus_name_loc(from_msg_key{1}),1);
        from_phases = bus_info(from).nodes';
        s_v = zeros(3,1);
        s_v(from_phases) = source_voltage;
        to_phases = bus_info(to).nodes';
        source_voltage = s_v(to_phases,1);
        ph_ff = zeros(3,1);
        ph_f = bus_name_loc(from_msg_key{1});
        ph_ff(from_phases) = ph_f;
        Wph = W(bus_name_loc(to_msg_key{1}),ph_ff(to_phases));
        all_voltage(bus_name_loc(to_msg_key{1}),1) = transpose(1/trace(source_voltage*source_voltage')*Wph*source_voltage);
        
    end
    
    figure(6)
    subplot(211)
    plot(abs(all_voltage(1:end)),'linewidth',1.5);
    hold on
    plot(abs(sorted_vcomplex(1:end)),'linewidth',1.5);
    diff = abs(sorted_vcomplex)-abs(all_voltage);
    legend('Estimated','Actual');
    fprintf('Difference in Voltage Magnitude Estimation:%f\n',max(abs(diff)));
    
    subplot(212)
    plot(angle(all_voltage(1:end)),'linewidth',1.5);
    
    hold on
    plot(angle(sorted_vcomplex(1:end)),'linewidth',1.5);
    diff = angle(sorted_vcomplex)-angle(all_voltage);
    legend('Estimated','Actual');
    fprintf('Difference in Voltage Angle Estimation:%f\n',max(abs(diff)));
    fprintf('Number of Sensors:%d\n',length(bus_name_sensor));
    
    
    [phase_a_buses_loc,~,phase_a_buses_values] = find(phase_distribution(:,1));
    [phase_b_buses_loc,~,phase_b_buses_values] = find(phase_distribution(:,2));
    [phase_c_buses_loc,~,phase_c_buses_values] = find(phase_distribution(:,3));
    
    abs_sorted_vcomplex = abs(sorted_vcomplex);
    abs_vcon = abs(all_voltage);
    phase_a_actual =  abs_sorted_vcomplex(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
    phase_a_vcon =  abs_vcon(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
    phase_b_actual =  abs_sorted_vcomplex(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
    phase_b_vcon =  abs_vcon(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
    phase_c_actual =  abs_sorted_vcomplex(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
    phase_c_vcon =  abs_vcon(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
    
    angle_sorted_vcomplex = angle(sorted_vcomplex);
    angle_vcon = angle(all_voltage);
    phase_a_actual_angle =  angle_sorted_vcomplex(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
    phase_a_vcon_angle =  angle_vcon(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
    phase_b_actual_angle =  angle_sorted_vcomplex(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
    phase_b_vcon_angle =  angle_vcon(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
    phase_c_actual_angle =  angle_sorted_vcomplex(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
    phase_c_vcon_angle =  angle_vcon(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
    
    
    
    figure(8)
    subplot(311)
    plot(phase_a_actual,'linewidth',1.5)
    hold on
    plot(phase_a_vcon,'linewidth',1.5)
    subplot(312)
    plot(phase_b_actual,'linewidth',1.5)
    hold on
    plot(phase_b_vcon,'linewidth',1.5)
    subplot(313)
    plot(phase_c_actual,'linewidth',1.5)
    hold on
    plot(phase_c_vcon,'linewidth',1.5)
    
    
    figure(9)
    subplot(311)
    plot(phase_a_actual_angle,'linewidth',1.5)
    hold on
    plot(phase_a_vcon_angle,'linewidth',1.5)
    subplot(312)
    plot(phase_b_actual_angle,'linewidth',1.5)
    hold on
    plot(phase_b_vcon_angle,'linewidth',1.5)
    subplot(313)
    plot(phase_c_actual_angle,'linewidth',1.5)
    hold on
    plot(phase_c_vcon_angle,'linewidth',1.5)
    
    
else
    fprintf('CVX_STATUES:%s\n',cvx_status);
    fprintf('Number of Sensors:%d\n',length(bus_name_sensor));
end

% csvwrite('phase_a_mag_low.csv',[transpose(1:length(phase_a_actual)) phase_a_actual phase_a_vcon])
% csvwrite('phase_b_mag_low.csv',[transpose(1:length(phase_b_actual)) phase_b_actual phase_b_vcon])
% csvwrite('phase_c_mag_low.csv',[transpose(1:length(phase_c_actual)) phase_c_actual phase_c_vcon])
% 
% 
% csvwrite('phase_a_angle_low.csv',[transpose(1:length(phase_a_actual_angle)) phase_a_actual_angle phase_a_vcon_angle])
% csvwrite('phase_b_angle_low.csv',[transpose(1:length(phase_b_actual_angle)) phase_b_actual_angle phase_b_vcon_angle])
% csvwrite('phase_c_angle_low.csv',[transpose(1:length(phase_c_actual_angle)) phase_c_actual_angle phase_c_vcon_angle])

