baseDir = '/media/adamb/bonczData/EEG_resting_state/alpha';
freq = 'alpha';
method = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
gammaParam = 1.05;
repNo = 1;
%groupNumbers = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
groupN = 10;
permNo = 1000;

% load data
connData = cell(length(method),1);
for m = 1:length(method)
    metric = method{m};
    connF = [baseDir, '/group_', freq, '_', metric, '.mat'];
    tmp = load(connF);
    connData{m} = tmp.connData;
end

% get sizes
[subNo, epochNo, roiNo, ~] = size(connData{1});

% preallocate for partitions
allParts = nan(permNo, length(method), roiNo);

% permutations
for permIdx = 1:permNo

    % define groups, random sampling
    groupSubs = randperm(subNo, groupN);

    for m = 1:length(method)
        
        % current data depends on metric
        data = connData{m};
        
        % subgroup average
        groupData = squeeze(mean(mean(data(groupSubs, :, :, :), 1), 2));

        % normalize
        groupData = normalizeMatrix(groupData, 'mean', false);

        % symm
        groupData = triu(groupData, 1) +  triu(groupData, 1)';

        % get B (null model)
        [B, ~] = modularity(groupData, gammaParam);

        % iterated louvain, repNo repetitions
        part = nan(repNo, roiNo);
        for n = 1:repNo
            part(n, :) = iterated_genlouvain(B, 5000, 0, [], 'moverandw');
        end

        if repNo > 1
            % get consensus partition based on nodal assoc matrix
            [S2, Q2, ~, ~] = consensus_iterative(part);
            % select best one
            tmpIdx = find(Q2==max(Q2), 1);
            % collect results
            group1part = S2(tmpIdx, :);
            allParts(permIdx, m, :) = S2(tmpIdx, :);

        else
            
            allParts(permIdx, m, :) = part;

        end
        
    end

end
    

% similarity

% preallocate for similarity
zrandRes = nan(permNo, length(method), length(method));

for p = 1:permNo
    permParts = squeeze(allParts(p, :, :));
    for m1 = 1:length(method)
        for m2 = 1:length(method)
            zrandRes(p, m1, m2) = zrand(permParts(m1, :), permParts(m2, :));
        end
    end
end

