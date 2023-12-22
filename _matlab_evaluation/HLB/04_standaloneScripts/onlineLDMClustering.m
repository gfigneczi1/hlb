function onlineLDMClustering(data, Pcentroids, clusters_in)

for i=1:length(data)
    P_est_Left = reshape(data(i).data(1).GT.Pest, 9,1);
    P_est_Right = reshape(data(i).data(2).GT.Pest, 9,1);
    P0 = data(i).p0;
    P_GT_matrix(i,:) = [P_est_Left' P_est_Right', P0];
end

%rng('default');
gt_clusters = kmeans(P_GT_matrix,3);
for i=1:3
    P_clusterizedMatrix = P_GT_matrix(gt_clusters==i,:);
    Pcentroids{i} = sum(P_clusterizedMatrix)/size(P_clusterizedMatrix,1);
end


    for i=1:length(data)
        P_est_Left = reshape(data(i).data(1).EKF.Pest, 9,1);
        P_est_Right = reshape(data(i).data(2).EKF.Pest, 9,1);
        P0 = data(i).p0;
        for j=1:length(Pcentroids)
            L2(i,j) = sum(sum(([P_est_Left; P_est_Right; P0]' - Pcentroids{j}).^2));
        end
        [~, clusters_out(i)] = min(L2(i,:));
    end
    fprintf("General clustering results:\n")
    disp(gt_clusters');
    fprintf("Estimated clustering results:\n")
    disp(clusters_out);
end

