function [sensor_p_all_included,sensor_p_all_excluded] = sensor_placement_matpower(mpc_case)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
define_constants;
f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

G=digraph;

for i = 1:length(f)
    G=addedge(G,f(i),t(i));
end


max_distance_node = 0;
max_distance = 0;
leaf_nodes = find(outdegree(G)==0);
sensor_p_all_included = leaf_nodes';
excluded_nodes = [];
for n=1:length(leaf_nodes)
    path = shortestpath(G,1,leaf_nodes(n));
    if length(path)>max_distance
        max_distance = length(path);
        max_distance_node = leaf_nodes(n);
    end
    for p = length(path):-2:1
        sensor_p_all_included = [sensor_p_all_included path(p)];
        if p >1
        excluded_nodes = [excluded_nodes path(p-1)];
        end
    end
    
end
sensor_p_all_included = unique(sensor_p_all_included);
zero_consumption_bus = find(mpc_case.bus(:,PD)==0);
sensor_p_all_included = setdiff(sensor_p_all_included,zero_consumption_bus);
sensor_p_all_included = [1 sensor_p_all_included];
sensor_p_all_included = unique(sensor_p_all_included);


excluded_nodes = unique(excluded_nodes);
sensor_p_all_excluded = 1:length(mpc_case.bus);
sensor_p_all_excluded(excluded_nodes) = [];

% max_distance_path = shortestpath(G,1,max_distance_node);
% sensor_p_all_excluded = unique([sensor_p_all_excluded max_distance_path]);

end

