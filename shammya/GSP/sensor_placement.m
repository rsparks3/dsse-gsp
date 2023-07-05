clc;
clear all;
close all;
define_constants;
casestr = 'case85.m';
mpc_case = loadcase(casestr);

f = mpc_case.branch(:,F_BUS);
t = mpc_case.branch(:,T_BUS);

G=digraph;

for i = 1:length(f)
    G=addedge(G,f(i),t(i));
end

hop_map = {};
%%
source = [1];
hop = 0;

while 1 
    new_source = [];
    for i = 1:length(source)
       new_source = [new_source  successors(G,source(i))'];
    end
    if isempty(new_source)
        break
    end
    hop = hop + 1 ;
    fprintf('%d-->',hop);
    disp(new_source);
    hop_map{hop} = new_source;
    fprintf('\n')
%     disp(new_source);
    
    source = new_source;
    
 
end

max_hop = hop;

sensor_locations = [1];

for hop = max_hop:-2:1
    sensor_locations = [sensor_locations hop_map{hop}];
end
sensor_locations = unique(sensor_locations);

% highlight(h,sensor_locations,'NodeColor','r')

%%


leaf_nodes = find(outdegree(G)==0);
sensor_p = leaf_nodes';
excluded_nodes = [];
for n=1:length(leaf_nodes)
    path = shortestpath(G,1,leaf_nodes(n));
    for p = length(path):-2:1
        sensor_p = [sensor_p path(p)];
        if p >1
        excluded_nodes = [excluded_nodes path(p-1)];
        end
    end
    
end

sensor_p = unique(sensor_p);
excluded_nodes = unique(excluded_nodes);
figure
h = plot(G);
highlight(h,sensor_p,'NodeColor','k')

figure
sensor_p_1 = 1:85;
sensor_p_1(excluded_nodes) = [];

h1 = plot(G);
highlight(h1,sensor_p_1,'NodeColor','k')




