load('D:\edgePruning_results\alpha\s02_alpha_edgePruningInfo.mat');

numberOfIterations = 100;

prunedConnectivity = prunedConn;
meanSurrogateData = meanSurrConn;

S_story1_array = [];
S_story2_array = [];
S_story3_array = [];
S_story4_array = [];
Q_story1_array = [];
Q_story2_array = [];
Q_story3_array = [];
Q_story4_array = [];
Q_story1_varianceArray = [];
Q_story2_varianceArray = [];
Q_story3_varianceArray = [];
Q_story4_varianceArray = [];
gamma_array = [];
omega_array = [];
tic
for omega = 0.0001 : 0.0001 % only one omega value (0.0001) is tested
    for gamma = 0.7 : 0.002 : 2.4
        [B] = calculateMultiLayerModularityMatrix(prunedConnectivity, meanSurrogateData, gamma, omega);
        Q_story1_max = 0;
        Q_story2_max = 0;
        Q_story3_max = 0;
        Q_story4_max = 0;
        Q_story1_arrayForVarianceCheck = zeros(1, numberOfIterations);
        Q_story2_arrayForVarianceCheck = zeros(1, numberOfIterations);
        Q_story3_arrayForVarianceCheck = zeros(1, numberOfIterations);
        Q_story4_arrayForVarianceCheck = zeros(1, numberOfIterations);
        for iterationIndex  = 1 : numberOfIterations
            [S_story1, Q_story1] = genlouvain(B{1});
            Q_story1_arrayForVarianceCheck(iterationIndex) = Q_story1;
            if Q_story1 > Q_story1_max
                Q_story1_max = Q_story1;
                S_story1_max = S_story1;
            end
            [S_story2, Q_story2] = genlouvain(B{2});
            Q_story2_arrayForVarianceCheck(iterationIndex) = Q_story2;
            if Q_story2 > Q_story2_max
                Q_story2_max = Q_story2;
                S_story2_max = S_story2;
            end
            [S_story3, Q_story3] = genlouvain(B{3});
            Q_story3_arrayForVarianceCheck(iterationIndex) = Q_story3;
            if Q_story3 > Q_story3_max
                Q_story3_max = Q_story3;
                S_story3_max = S_story3;
            end
            [S_story4, Q_story4] = genlouvain(B{4});
            Q_story4_arrayForVarianceCheck(iterationIndex) = Q_story4;
            if Q_story4 > Q_story4_max
                Q_story4_max = Q_story4;
                S_story4_max = S_story4;
            end
        end
        S_story1_array = [S_story1_array S_story1_max];
        S_story2_array = [S_story2_array S_story2_max];
        S_story3_array = [S_story3_array S_story3_max];
        S_story4_array = [S_story4_array S_story4_max];
        Q_story1_array = [Q_story1_array Q_story1_max];
        Q_story2_array = [Q_story2_array Q_story2_max];
        Q_story3_array = [Q_story3_array Q_story3_max];
        Q_story4_array = [Q_story4_array Q_story4_max];
        Q_story1_varianceArray = [Q_story1_varianceArray var(Q_story1_arrayForVarianceCheck)];
        Q_story2_varianceArray = [Q_story2_varianceArray var(Q_story2_arrayForVarianceCheck)];
        Q_story3_varianceArray = [Q_story3_varianceArray var(Q_story3_arrayForVarianceCheck)];
        Q_story4_varianceArray = [Q_story4_varianceArray var(Q_story4_arrayForVarianceCheck)];
        gamma_array = [gamma_array gamma];
        omega_array = [omega_array omega];
    end
end
toc
numberOfCommunities_story1 = max(S_story1_array);
numberOfCommunities_story2 = max(S_story2_array);
numberOfCommunities_story3 = max(S_story3_array);
numberOfCommunities_story4 = max(S_story4_array);

figure();
plot(gamma_array, numberOfCommunities_story1);
xlabel('Structural resolution parameter - gamma');
ylabel('Number of communities');
title('Story 1');
figure();
plot(gamma_array, numberOfCommunities_story2);
xlabel('Structural resolution parameter - gamma');
ylabel('Number of communities');
title('Story 2');
figure();
plot(gamma_array, numberOfCommunities_story3);
xlabel('Structural resolution parameter - gamma');
ylabel('Number of communities');
title('Story 3');
figure();
plot(gamma_array, numberOfCommunities_story4);
xlabel('Structural resolution parameter - gamma');
ylabel('Number of communities');
title('Story 4');

figure();
plot(gamma_array, Q_story1_varianceArray);
xlabel('Structural resolution parameter - gamma');
ylabel('Variance of Q');
title('Story 1');
figure();
plot(gamma_array, Q_story2_varianceArray);
xlabel('Structural resolution parameter - gamma');
ylabel('Variance of Q');
title('Story 2');
figure();
plot(gamma_array, Q_story3_varianceArray);
xlabel('Structural resolution parameter - gamma');
ylabel('Variance of Q');
title('Story 3');
figure();
plot(gamma_array, Q_story4_varianceArray);
xlabel('Structural resolution parameter - gamma');
ylabel('Variance of Q');
title('Story 4');

JaccardDistances_story1 = [];
for gammaValueIndex = 1 : size(S_story1_array, 2)-1
    type1_counter = 0;
    type2_counter = 0;
    type3_counter = 0;
    partitionVector1 = S_story1_array(:, gammaValueIndex);
    partitionVector2 = S_story1_array(:, gammaValueIndex+1);
    for firstNodeIndex = 1 : numel(partitionVector1)-1
        for secondNodeIndex = firstNodeIndex+1 : numel(partitionVector1)
            if (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type1_counter = type1_counter + 1;
            elseif (partitionVector1(firstNodeIndex) ~= partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type2_counter = type2_counter + 1;
            elseif (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) ~= partitionVector2(secondNodeIndex))
                type3_counter = type3_counter + 1;
            end
        end
    end
    JaccardIndex = type1_counter / (type1_counter + type2_counter + type3_counter);
    JaccardDistances_story1 = [JaccardDistances_story1 (1-JaccardIndex)];
end

JaccardDistances_story2 = [];
for gammaValueIndex = 1 : size(S_story2_array, 2)-1
    type1_counter = 0;
    type2_counter = 0;
    type3_counter = 0;
    partitionVector1 = S_story2_array(:, gammaValueIndex);
    partitionVector2 = S_story2_array(:, gammaValueIndex+1);
    for firstNodeIndex = 1 : numel(partitionVector1)-1
        for secondNodeIndex = firstNodeIndex+1 : numel(partitionVector1)
            if (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type1_counter = type1_counter + 1;
            elseif (partitionVector1(firstNodeIndex) ~= partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type2_counter = type2_counter + 1;
            elseif (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) ~= partitionVector2(secondNodeIndex))
                type3_counter = type3_counter + 1;
            end
        end
    end
    JaccardIndex = type1_counter / (type1_counter + type2_counter + type3_counter);
    JaccardDistances_story2 = [JaccardDistances_story2 (1-JaccardIndex)];
end

JaccardDistances_story3 = [];
for gammaValueIndex = 1 : size(S_story3_array, 2)-1
    type1_counter = 0;
    type2_counter = 0;
    type3_counter = 0;
    partitionVector1 = S_story3_array(:, gammaValueIndex);
    partitionVector2 = S_story3_array(:, gammaValueIndex+1);
    for firstNodeIndex = 1 : numel(partitionVector1)-1
        for secondNodeIndex = firstNodeIndex+1 : numel(partitionVector1)
            if (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type1_counter = type1_counter + 1;
            elseif (partitionVector1(firstNodeIndex) ~= partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type2_counter = type2_counter + 1;
            elseif (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) ~= partitionVector2(secondNodeIndex))
                type3_counter = type3_counter + 1;
            end
        end
    end
    JaccardIndex = type1_counter / (type1_counter + type2_counter + type3_counter);
    JaccardDistances_story3 = [JaccardDistances_story3 (1-JaccardIndex)];
end

JaccardDistances_story4 = [];
for gammaValueIndex = 1 : size(S_story4_array, 2)-1
    type1_counter = 0;
    type2_counter = 0;
    type3_counter = 0;
    partitionVector1 = S_story4_array(:, gammaValueIndex);
    partitionVector2 = S_story4_array(:, gammaValueIndex+1);
    for firstNodeIndex = 1 : numel(partitionVector1)-1
        for secondNodeIndex = firstNodeIndex+1 : numel(partitionVector1)
            if (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type1_counter = type1_counter + 1;
            elseif (partitionVector1(firstNodeIndex) ~= partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) == partitionVector2(secondNodeIndex))
                type2_counter = type2_counter + 1;
            elseif (partitionVector1(firstNodeIndex) == partitionVector1(secondNodeIndex)) && (partitionVector2(firstNodeIndex) ~= partitionVector2(secondNodeIndex))
                type3_counter = type3_counter + 1;
            end
        end
    end
    JaccardIndex = type1_counter / (type1_counter + type2_counter + type3_counter);
    JaccardDistances_story4 = [JaccardDistances_story4 (1-JaccardIndex)];
end

figure();
plot(gamma_array(1:end-1), JaccardDistances_story1);
xlabel('Structural resolution parameter - gamma');
ylabel('Jaccard distance');
title('Story 1');
figure();
plot(gamma_array(1:end-1), JaccardDistances_story2);
xlabel('Structural resolution parameter - gamma');
ylabel('Jaccard distance');
title('Story 2');
figure();
plot(gamma_array(1:end-1), JaccardDistances_story3);
xlabel('Structural resolution parameter - gamma');
ylabel('Jaccard distance');
title('Story 3');
figure();
plot(gamma_array(1:end-1), JaccardDistances_story4);
xlabel('Structural resolution parameter - gamma');
ylabel('Jaccard distance');
title('Story 4');

save('community_workspace_variables.mat');

% N = 62;
% T = 40;
% S1_max = reshape(S1_max,N,T);
% S2_max = reshape(S1_max,N,T);
% S3_max = reshape(S1_max,N,T);
% S4_max = reshape(S1_max,N,T);