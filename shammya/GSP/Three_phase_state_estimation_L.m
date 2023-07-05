clear all;
clc;

load('123bus.mat')

%%
Uk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
YUk_UkH =  zeros(number_of_eigenvectors,number_of_eigenvectors*number_of_eigenvectors);
A=U_k;
C=U_k';
for k = 1:nb
   Uk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end
A=Y_DD*U_k;
for k = 1:nb
   YUk_UkH(k,:) =  transpose(vec(C(:,k)*A(k,:)));
end

cvx_clear
cvx_solver SDPT3
cvx_precision low
cvx_begin 

    variable tildeW(number_of_eigenvectors,number_of_eigenvectors) complex semidefinite
    variable s_est(nb) complex
    variable v_est(nb)
    variable a_est(nb)
    minimize (15*norm(sm - M*s_est, 2) + 5*norm(vm.*vm-M*(v_est),2) + 5*norm(cm.*cm-M*a_est,2) )

    subject to
    %         sm == M*s_est;
    tildeW == tildeW';
%     s_est == conj(diag(Y_DD*U_k*tildeW*U_k'));
%     v_est == diag(U_k*tildeW*(U_k)');
    for k = 1:nb
        v_est(k) == Uk_UkH(k,:)*vec(transpose(tildeW)); 
        s_est(k) == conj(YUk_UkH(k,:)*vec(transpose(tildeW))); 
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
    diff = abs(sorted_vcomplex)-abs(all_voltage);
    fprintf('Difference in Voltage Magnitude Estimation:%f\n',max(abs(diff)));
    diff = angle(sorted_vcomplex)-angle(all_voltage);
    fprintf('Difference in Voltage Angle Estimation:%f\n',max(abs(diff)));
    fprintf('Number of Sensors:%d\n',length(bus_name_sensor));
%     figure(6)
%     subplot(211)
%     plot(abs(all_voltage(1:end)),'linewidth',1.5);
%     hold on 
%     plot(abs(sorted_vcomplex(1:end)),'linewidth',1.5);
%     legend('Estimated','Actual');

% 
%     subplot(212)
%     plot(angle(all_voltage(1:end)),'linewidth',1.5);
%     hold on 
%     plot(angle(sorted_vcomplex(1:end)),'linewidth',1.5);
%     
%     legend('Estimated','Actual');
%     
%     
%     [phase_a_buses_loc,~,phase_a_buses_values] = find(phase_distribution(:,1));
%     [phase_b_buses_loc,~,phase_b_buses_values] = find(phase_distribution(:,2));
%     [phase_c_buses_loc,~,phase_c_buses_values] = find(phase_distribution(:,3));
%     
%     abs_sorted_vcomplex = abs(sorted_vcomplex);
%     abs_vcon = abs(all_voltage);
%     phase_a_actual =  abs_sorted_vcomplex(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
%     phase_a_vcon =  abs_vcon(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
%     phase_b_actual =  abs_sorted_vcomplex(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
%     phase_b_vcon =  abs_vcon(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
%     phase_c_actual =  abs_sorted_vcomplex(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
%     phase_c_vcon =  abs_vcon(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
%     
%     angle_sorted_vcomplex = angle(sorted_vcomplex);
%     angle_vcon = angle(all_voltage);
%     phase_a_actual_angle =  angle_sorted_vcomplex(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
%     phase_a_vcon_angle =  angle_vcon(phase_a_buses_values)./base_voltage(phase_a_buses_loc) ;
%     phase_b_actual_angle =  angle_sorted_vcomplex(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
%     phase_b_vcon_angle =  angle_vcon(phase_b_buses_values)./base_voltage(phase_b_buses_loc) ;
%     phase_c_actual_angle =  angle_sorted_vcomplex(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
%     phase_c_vcon_angle =  angle_vcon(phase_c_buses_values)./base_voltage(phase_c_buses_loc) ;
%     
%     
%     
%     figure(8)
%     subplot(311)
%     plot(phase_a_actual,'linewidth',1.5)
%     hold on 
%     plot(phase_a_vcon,'linewidth',1.5)
%     subplot(312)
%     plot(phase_b_actual,'linewidth',1.5)
%     hold on 
%     plot(phase_b_vcon,'linewidth',1.5)
%     subplot(313)
%     plot(phase_c_actual,'linewidth',1.5)
%     hold on 
%     plot(phase_c_vcon,'linewidth',1.5)
%     
%     
%     figure(9)
%     subplot(311)
%     plot(phase_a_actual_angle,'linewidth',1.5)
%     hold on 
%     plot(phase_a_vcon_angle,'linewidth',1.5)
%     subplot(312)
%     plot(phase_b_actual_angle,'linewidth',1.5)
%     hold on 
%     plot(phase_b_vcon_angle,'linewidth',1.5)
%     subplot(313)
%     plot(phase_c_actual_angle,'linewidth',1.5)
%     hold on 
%     plot(phase_c_vcon_angle,'linewidth',1.5)
else
    fprintf('CVX_STATUES:%s\n',cvx_status);
    fprintf('Number of Sensors:%d\n',length(bus_name_sensor));
end
