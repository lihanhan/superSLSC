function alllabels_result = SVLA_supervoxel(allpoints,lamda1,lamda2,resolution)
    
    alllabels = allpoints(:,end);
    %% 构建KDtree，得到每个超体素的邻接超体素
    allpsKdtree = KDTreeSearcher(allpoints(:,1:3));
    nearpoints_num = 5; % 选取邻域点的数量
    allfd=[];   % 存储每个点粗糙性的容器
    all_normal=[];  % 存储每个点法向量的容器
    all_neipsid=[];     % 存储每个点邻域的容器
    
    for i = 1:size(allpoints,1)
        Idx = knnsearch(allpsKdtree, allpoints(i,1:3),'k',10);% 搜索最近点较慢
        [value, normal_vector] = pca_function(allpoints(Idx,1:3));
        all_neipsid = [all_neipsid;Idx];
        all_normal = [all_normal;normal_vector];
        tempfd = sqrt(value(3))/sqrt(value(1));
        allfd = [allfd; tempfd];
    end
    %% 计算每个超体素块所属的极值点
    un_labels = unique(alllabels,'stable');
    % 极值点的id
    all_extreme_ps = [];
    for i=1:size(un_labels,1)
        pointids = find(alllabels == un_labels(i));% 超体素内点的id
        t_super_fd = allfd(pointids);% 超体素内点的fd
        % 当前超体素的极值点
        extreme_pid = pointids(find(min(t_super_fd)==t_super_fd));
        if length(extreme_pid)>1
            extreme_pid = extreme_pid(1);
        end
        all_extreme_ps = [all_extreme_ps;extreme_pid];
    end
    ExtremeKdtree = KDTreeSearcher(allpoints(all_extreme_ps,1:3));% 极值点树
    
    %% nei_input为极值点xyz-label
    nei_input = [allpoints(all_extreme_ps,1:3),alllabels(all_extreme_ps)];
    %% 得到每个超体素的邻接超体素,保证每个输入点均有邻域点信息
    %% 自适应地求邻域点信息
%    Glabel_result = cell([],1);
%    fresult = Judge_empty(Glabel_result);
%     while fresult
%         [Glabel_result,DT_index ,DT_label] = Get_Adjacent_TIN(nei_input,alpha_value);
%         alpha_value = alpha_value + 1;
%         fresult = Judge_empty(Glabel_result);
%     end
    %% 求每个超体素的邻域label
    exIdxlabel =[];
    for i=1:size(nei_input)
        texIdx = knnsearch(ExtremeKdtree, nei_input(i,1:3),'k',4);
        t_nei=[];
        for j = 1:length(texIdx)
            t_nei = [t_nei,nei_input(texIdx(j),end)];
        end
        exIdxlabel =[exIdxlabel;t_nei];
    end
    
    
    Travel = zeros(length(allpoints),1);
    for i=1:size(allpoints,1)
        if Travel(i)==1
            continue;
        end
        t_pointlabel= alllabels(i);   % 当前点所属的超体素label
        % 得到i点所属超体素的邻域
        t_p_l_id = find(t_pointlabel==un_labels(:));
        %nei_ids = Glabel_result{t_p_l_id};
        %% 得到邻域点的label
        %% 计算Fc
        t_allnei_ex_ps=[];  % 与第i个点相邻的 所有超体素极值点
        exIdxlabel_=exIdxlabel(t_p_l_id,:);
        for j=1:length(exIdxlabel_)
            t_neilabel = exIdxlabel_(j);
            t_neilabel_id = find(t_neilabel==un_labels);
            t_super_ex_ps = all_extreme_ps(t_neilabel_id,:);
            t_allnei_ex_ps=[t_allnei_ex_ps;t_super_ex_ps];
        end
        % 得到最近的极值点id
        resultid = Cal_Ci(allpoints(i,1:3),t_allnei_ex_ps,allpsKdtree,allpoints,lamda1,lamda2,resolution);
        % 当前超体素内的点
        tempsuper_id = find(alllabels==t_pointlabel);
        Travel(tempsuper_id) = 1;
        % 将第i个点的label 变为  最近极值点的label
        alllabels(tempsuper_id) = alllabels(resultid);
        %alllabels(tempsuper_id) = un_labels(resultid); %%---------此处是否是un_labels
    end
    alllabels_result = alllabels;
    
end