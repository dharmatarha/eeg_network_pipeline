% quick and dirty script to load the resting state data of the
% listener-listener hyperscan subjects

% Basic params needed
connMethod = 'iplv';
freq = 'alpha';
subs = {'s02','s03','s04','s05','s06','s07','s08','s09'...
           ,'s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
           's21','s22','s23','s24','s25','s26','s27','s28'};

% loaction of subject-level resting-state connectivity files
connFold = ['/media/adamb/bonczData/EEG_resting_state/', freq, '/'];

% load already averaged resting state data 
connData = load([connFold, 'group_', freq, '_', connMethod, '.mat']);

% only subject s22 is missing from the RS data, just grab the first 28
% subjects' data
RS_subs = connData.subjects(1:28);
connData = connData.connData(1:28, :, :, :);

% average across epochs for subjects of interest
meanRS = squeeze(mean(connData, 2));

% delete data from l3_s01, l3_s10 and l3_s29 (they are not part of final
% hyperscan data set)
delSubs = {'l3_s01', 'l3_s10', 'l3_s29'};
for i = 1:3
    delS = delSubs{i};
    delMask = ismember(RS_subs, delS);
    meanRS(delMask, :, :) = [];
    RS_subs(delMask) = [];
end
    
% group average RS
groupRS = squeeze(mean(meanRS, 1));


% % get file path for each file we are interested in
% connFiles = cell(numel(subs), 1);
% for s = 1:numel(subs)
%     subID = subs{s};
%     connFiles{s} = [connFold, 'l3_', subID, '_', freq, '_', connMethod, '.mat'];
%     if ~exist(connFiles{s})
%         error('Cannot find subject-level connectivity file!');
%     end
% end


       
       