function final_labels = Exchange(points,neinum)

    labels = points(:,end);
    firstlabel = unique(points(:,end),'stable');
    %% 计算每个label超体素的代表中心点id
    represent_id = [];
    for i=1:size(firstlabel,1)
        tlabel = firstlabel(i);
        tempsuperpoints_id = find(tlabel==points(:,end));
        % 计算当前超体素的中心点坐标
        temppoint = Calculate_center(points(tempsuperpoints_id,1:3));   % 以距离为划分，得到的中心
        tKdtree = KDTreeSearcher(points(tempsuperpoints_id,1:3));
        tIdx = knnsearch(tKdtree, temppoint,'k',1);
        tcenterIdx = tempsuperpoints_id(tIdx);
        represent_id = [represent_id; tcenterIdx];
    end
    %% 计算每个点到其所属中心点的距离
    dis=[];
    for i=1:size(points,1)
        t_l = points(i,end);
        t_lid = find(firstlabel==t_l);
        t_re_pid = represent_id(t_lid);% 当前label对应的种子点的id
        if i~=t_re_pid
            t_d = Cal_sim_dis(points,i,t_re_pid);
        else
            t_d = 0;
        end
        dis = [dis;t_d];
    end
    
    %% 获取每个点的邻域
    allpsKdtree = KDTreeSearcher(points(:,1:3));
    allpoints_nei=[];
    
    queue=[];
    in_q = zeros(size(points,1),1);
    %all_normals = [];
    for i=1:length(points)
        t_neiIdx = knnsearch(allpsKdtree, points(i,1:3),'k',neinum);
        %t_normal = pca(points(t_neiIdx,1:3));
        %all_normals =[all_normals;t_normal];
        for j=1:size(t_neiIdx,2)
            t_p_id = t_neiIdx(j);
            if points(t_p_id,end)~=points(i,end)
                if in_q(i)==0
                    queue = [queue;i];
                    in_q(i)=1;
                end
                if in_q(j)==0
                    queue = [queue;j];
                    in_q(j)=1;
                end
            end
        end
        allpoints_nei=[allpoints_nei;t_neiIdx];
    end
    %% 开始交换
    while(~isempty(queue))
        % 提取队列的第一个点
        t_id = queue(1); queue(1) = [];
        in_q(t_id) = 0;
        change = 0;
        for i=1:length(allpoints_nei(t_id,:))
            a = points(t_id,end);% 当前点的label
            t_j = allpoints_nei(t_id,i);
            b = points(t_j,end);% 邻域点的label
            if a==b
                continue;
            end
            % 计算第i个点到第b个label代表点的 距离
            t_b_id = find(firstlabel==b);
            nei_id = represent_id(t_b_id);
            d = Cal_sim_dis(points,t_id,nei_id);     % represent_id(t_b_id)为代表点的id
            if (d < dis(t_id)) 
                labels(t_id) = b;
                dis(t_id) = d;
                change = 1;
            end
        end
        if change
            % 将邻域点放置在队列中,继续交换
            for j=1:size(allpoints_nei(t_id,:))
                t_j = allpoints_nei(t_id,i);
                if  labels(t_id) ~= labels(t_j)
                    if ~in_q(t_j)
                        queue = [queue; t_j];
                        in_q(t_j) = 1;
                    end
                end
            end
        end
    end
    final_labels = labels;
end