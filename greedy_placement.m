function [P_S,M_s,ammeter_lines] = greedy_placement(U_K,M,M_s,location_choice,G,f,t,method,Yuk,Yik)
%% greedy implementation of optimal placement:
% M_s is a set of already present indices/locations-just need to add onto it
% each time-column vector
% [f,t] = findedge(G);
if method == 1
    M_tilde = M-length(M_s); % total number to be found
    [N,~]=size(U_K);
    i=1;
    
    node_degress = outdegree(G,location_choice);
    child_node_loc = node_degress==0;
    location_choice(child_node_loc) = [];
    ammeter_lines = [];
    j_tilde = location_choice;
    j_tilde(M_s)=[];
    while i <=M_tilde
        %sigma_min=0;sigma_ind=1;
        sigma_min=zeros(length(j_tilde),1);
        for j=1: length(j_tilde)
            A = U_K( [M_s;j_tilde(j)] ,:);
            [~,sigma_min(j),~,flag] = svds(A,1,'smallest');
            if(flag)
                continue;
            end
            % sigma_temp = svds(A,1,‘largest’); %min the largest
            %       if(sigma_temp>sigma_min)
            %       %if(sigma_temp<sigma_min)
            %         sigma_min = sigma_temp;
            %         sigma_ind=j;
            %       end
        end
        [~,min_ind]=max(sigma_min);
        M_s = [M_s;j_tilde(min_ind)];
        j_tilde(min_ind)=[];
        i=i+1;
    end
    P_S = eye(N);
    P_S=P_S(:,M_s');
    
    for i = 1:length(M_s)
        node = M_s(i);
        from_node_loc = f==node;
        to_nodes = t (from_node_loc);
        node_degress = outdegree(G,to_nodes);
        [~,index_max] = max(node_degress);
        ammeter_lines = [ammeter_lines find(t==to_nodes(index_max))];
        
    end
end

if method == 2
    M_tilde = M-length(M_s); % total number to be found
    [N,~]=size(U_K);
    i=1;
    ammeter_lines = [];
    j_tilde = location_choice;
    j_tilde(M_s)=[];
    while i <=M_tilde
        %sigma_min=0;sigma_ind=1;
        sigma_min=zeros(length(j_tilde),1);
        for j=1: length(j_tilde)
            A = U_K([M_s;j_tilde(j)],:);
            [~,sigma_min(j),~,flag] = svds(A,1,'smallest');
            if(flag)
                continue;
            end
            % sigma_temp = svds(A,1,‘largest’); %min the largest
            %       if(sigma_temp>sigma_min)
            %       %if(sigma_temp<sigma_min)
            %         sigma_min = sigma_temp;
            %         sigma_ind=j;
            %       end
        end
        [~,min_ind]=max(sigma_min);
        M_s = [M_s;j_tilde(min_ind)];
        j_tilde(min_ind)=[];
        i=i+1;
    end
    P_S = eye(N);
    P_S=P_S(:,M_s');
    
    for i = 1:length(M_s)
        node = M_s(i);
        if outdegree(G,node) == 0
            ammeter_lines = [ammeter_lines find(t==node)];
            
        else
            from_node_loc = f==node;
            to_nodes = t (from_node_loc);
            node_degress = outdegree(G,to_nodes);
            [~,index_max] = max(node_degress);
            ammeter_lines = [ammeter_lines find(t==to_nodes(index_max))];
        end
        
    end
end


if method == 3
    M_tilde = M-length(M_s); % total number to be found
    [N,~]=size(U_K);
    i=1;
    node_degress = outdegree(G,location_choice);
    child_node_loc = node_degress==0;
    location_choice(child_node_loc) = [];
    ammeter_lines = [];
    j_tilde = location_choice;
    j_tilde(M_s)=[];
    while i <= M_tilde
        sigma_min_uk=zeros(length(j_tilde),1);
        sigma_min_yuk=zeros(length(j_tilde),1);
        for j=1: length(j_tilde)
            A = U_K([M_s;j_tilde(j)],:);
            [~,sigma_min_uk(j),~,flag] = svds(A,1,'smallest');
            if(flag)
                continue;
            end
            B = Yuk([M_s;j_tilde(j)],:);
            [~,sigma_min_yuk(j),~,flag] = svds(B,1,'smallest');
            if(flag)
                continue;
            end
        end
%         [~,min_ind]=max(abs(sigma_min_uk)+abs(sigma_min_yuk));
        [~,min_ind]= max(min([sigma_min_uk sigma_min_yuk],[],2));
        M_s = [M_s;j_tilde(min_ind)];
        j_tilde(min_ind)=[];
        i=i+1;
    end
    P_S = eye(N);
    P_S=P_S(:,M_s');
    
    for i = 1:length(M_s)
        node = M_s(i);
        from_node_loc = f==node;
        to_nodes = t (from_node_loc);
        node_degress = outdegree(G,to_nodes);
        [~,index_max] = max(node_degress);
        ammeter_lines = [ammeter_lines find(t==to_nodes(index_max))];
        
    end
  
end


if method == 4

    M_tilde = M-length(M_s); % total number to be found
    [N,~]=size(U_K);
    i=1;
%     node_degress = outdegree(G,location_choice);
%     child_node_loc = node_degress==0;
%     location_choice(child_node_loc) = [];
    ammeter_lines = [];
    j_tilde = location_choice;
    j_tilde(M_s)=[];
    while i <= M_tilde
        sigma_min_uk=zeros(length(j_tilde),1);
        sigma_min_yuk=zeros(length(j_tilde),1);
        for j=1: length(j_tilde)
            A = U_K([M_s;j_tilde(j)],:);
            [~,sigma_min_uk(j),~,flag] = svds(A,1,'smallest');
            if(flag)
                continue;
            end
            B = Yuk([M_s;j_tilde(j)],:);
            [~,sigma_min_yuk(j),~,flag] = svds(B,1,'smallest');
            if(flag)
                continue;
            end
        end
%         [~,min_ind]=max(abs(sigma_min_uk)+abs(sigma_min_yuk));
        [~,min_ind]= max(min([sigma_min_uk sigma_min_yuk],[],2));
        
        M_s = [M_s;j_tilde(min_ind)];
        j_tilde(min_ind)=[];
        i=i+1;
    end
    P_S = eye(N);
    P_S=P_S(:,M_s');
    
 
end





end