clc;
clear all;
% close all;


no_cluster = load('case85_no_clustering_scenario.mat');
% no_cluster = load('case85_clustering_scenario_noloss.mat');
actual_angle = no_cluster.actual_angle;
actual_voltage = no_cluster.actual_voltage;
estimated_angle = no_cluster.estimated_angle;
estimated_voltage = no_cluster.estimated_voltage;
failed_scenario = no_cluster.failed_scenario;

actual_angle(:,failed_scenario) = [];
actual_voltage(:,failed_scenario) = [];
estimated_angle(:,failed_scenario) = [];
estimated_voltage(:,failed_scenario) = [];

rms_angle_error = rms(estimated_angle-actual_angle);
[~,index] = sort(rms_angle_error,'ascend');
index = index(1:5000);

angle_error_plot = abs(actual_angle(:,index)-estimated_angle(:,index));
magnitude_error_plot = abs(actual_voltage(:,index)-estimated_voltage(:,index));

figure(5)
set(gca,'FontSize',18)

subplot(211)

imagesc(magnitude_error_plot)
title('No Clustering Only AMI','FontSize', 24)
c = colorbar;
c.FontSize=20;
ylabel('Magnitude Error (pu)','FontSize', 24)

subplot(212)
imagesc(angle_error_plot)
title('No Clustering Only AMI','FontSize', 24)
c = colorbar;
c.FontSize=20;
ylabel('Angle Error Plot','FontSize', 24)

%%
casestr = 'case85';
opt = mpoption('verbose',0,'out.all',0);
mpc_case = loadcase(casestr); % load case to get back to the original network
define_constants;
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);
edge_weights =  sqrt(mpc_case.branch(:,BR_R).^2 +  mpc_case.branch(:,BR_X).^2);
G= graph(f,t,edge_weights);

mag_error = abs(actual_voltage(:,index)-estimated_voltage(:,index));
avg_mag_error = mean(mag_error,2);

angle_error = abs(actual_angle(:,index)-estimated_angle(:,index));
avg_angle_error = mean(angle_error,2);

sensor_locations = load('case85_clustering_scenario_noloss.mat');
sensor_locations = sensor_locations.sensor_locations;
close all
%%
figure(5)
subplot(121)
p = plot(G);
cmap = jet;
colormap(gca, cmap);
G.Nodes.value =  avg_mag_error;
G.Nodes.NodeColors = G.Nodes.value;
p.NodeCData = G.Nodes.NodeColors;
Cdata = p.NodeCData;
 cmin = min(Cdata(:));
 cmax = max(Cdata(:));
 m = length(cmap);
 index = fix((Cdata-cmin)/(cmax-cmin)*m)+1;
q= plot(G,'NodeCData',G.Nodes.NodeColors);
labelnode(q,sensor_locations,'s')
labelnode(q,setdiff(1:85,sensor_locations),'')
RGB = squeeze(ind2rgb(index,cmap));
q.NodeFontSize = 14;
title('Avg Magnitude Error')
c = colorbar;
% colorbar(cmap)
fileID = fopen('vm.txt','w');
for i = 1:85
   data = strcat('\\definecolor{c',string(i),'}{rgb}{',string(RGB(i,1)),',',string(RGB(i,2)),',',string(RGB(i,3)),'}');
   fprintf(fileID,data) ;
   fprintf(fileID,'\n');
end
fclose(fileID);

subplot(122)
p = plot(G);
G.Nodes.value =  avg_angle_error;
G.Nodes.NodeColors = G.Nodes.value;
p.NodeCData = G.Nodes.NodeColors;
Cdata = p.NodeCData;
 cmin = min(Cdata(:));
 cmax = max(Cdata(:));
 m = length(cmap);
 index = fix((Cdata-cmin)/(cmax-cmin)*m)+1;
q= plot(G,'NodeCData',G.Nodes.NodeColors);
labelnode(q,sensor_locations,'s')
labelnode(q,setdiff(1:85,sensor_locations),'')
RGBa = squeeze(ind2rgb(index,cmap));
q.NodeFontSize = 14;
title('Avg Angle Error')
colorbar

fileID = fopen('va.txt','w');
for i = 1:85
   data = strcat('\\definecolor{c',string(i),'}{rgb}{',string(RGBa(i,1)),',',string(RGBa(i,2)),',',string(RGBa(i,3)),'}');
   fprintf(fileID,data) ;
   fprintf(fileID,'\n');
end
fclose(fileID);
