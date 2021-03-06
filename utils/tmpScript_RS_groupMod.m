baseDir = '/media/adamb/bonczData/EEG_resting_state/alpha';
freq = 'alpha';
method = 'orthAmpCorr';
gammaParam = 1.00;
repNo = 1;
groupNumbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
permNo = 1000;

% load data
connF = [baseDir, '/group_', freq, '_', method, '.mat'];
tmp = load(connF);
connData = tmp.connData;

% get sizes
[subNo, epochNo, roiNo, ~] = size(connData);

% preallocate for similarity
zrandRes = nan(permNo, length(groupNumbers));

for groupIdx = 1:length(groupNumbers)
    groupN = groupNumbers(groupIdx);

    % user message
    disp(['working on groupN = ', num2str(groupN)]);
    
    % permutations
    for permIdx = 1:permNo

        % define groups, random sampling
        groupSubs = randperm(subNo, groupN*2);
        group1 = groupSubs(1:groupN);
        group2 = groupSubs(groupN+1:end);

        % subgroup averages
        groupData1 = squeeze(mean(mean(connData(group1, :, :, :), 1), 2));
        groupData2 = squeeze(mean(mean(connData(group2, :, :, :), 1), 2));

        % normalize
        groupData1 = normalizeMatrix(groupData1, 'mean', false);
        groupData2 = normalizeMatrix(groupData2, 'mean', false);

        % symm
        groupData1 = triu(groupData1, 1) +  triu(groupData1, 1)';
        groupData2 = triu(groupData2, 1) +  triu(groupData2, 1)';

        % get B (null model)
        [B1, ~] = modularity(groupData1, gammaParam);
        [B2, ~] = modularity(groupData2, gammaParam);

        % iterated louvain, repNo repetitions
        part1 = nan(repNo, roiNo);
        part2 = nan(repNo, roiNo);
        for n = 1:repNo
            part1(n, :) = iterated_genlouvain(B1, 5000, 0, [], 'moverandw');
            part2(n, :) = iterated_genlouvain(B2, 5000, 0, [], 'moverandw');
        end

        if repNo > 1
            % get consensus partition based on nodal assoc matrix
            [S2, Q2, ~, ~] = consensus_iterative(part1);
            % select best one
            tmpIdx = find(Q2==max(Q2), 1);
            % collect results
            group1part = S2(tmpIdx, :);

            % get consensus partition based on nodal assoc matrix
            [S2, Q2, ~, ~] = consensus_iterative(part2);
            % select best one
            tmpIdx = find(Q2==max(Q2), 1);
            % collect results
            group2part = S2(tmpIdx, :);

            zrandRes(permIdx, groupIdx) = zrand(group1part, group2part);

        else
            zrandRes(permIdx, groupIdx) = zrand(part1, part2);

        end

    end
    
end




