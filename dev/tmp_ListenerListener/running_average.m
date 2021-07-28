function running_mean_upd = running_average(running_mean, seq_n, new_data)
%% Running average calculator
%
% USAGE: running_mean_upd = runnning_average(running_mean, seq_n, new_data)
%
% Set the first sample / dataset as the initial running average, then call
% it repeatedly with each new sample / dataset. 
% NOTE: After setting the
% first sample as running average, "seq_n" must equal 2 when calling
% "running_average" the first time, with the second sample as "new_data".
%
% Inputs:
% running_mean      - Numeric array. Running mean so far. 
% seq_n             - Numeric value, positive integer, sequential number 
%                   of the data to be averaged. 
%                   E.g. if in the current run you add the 10th
%                   sample to the running average, seq_n = 10.
% new_data          - Numeric array, data to be averaged. Same size as
%                   "running_mean".
%
% Outputs:
% running_mean_upd  - Numeric array, the running average updated with
%                   "new_data".
%


%% Input checks

if nargin ~= 3
    error(['Function "running_average" requires inputs "running_mean", ',...
        '"seq_n" and "new_data"!']);
end
if ~isnumeric(running_mean)
    error('Input arg "running_mean" should be a numeric matrix!');
end
if ~isnumeric(new_data) || ~isequal(size(running_mean), size(new_data))
    error('Input arg "new_data" should be a numeric matrix with the same size as "running_mean"!');
end
if ~isnumeric(seq_n) || numel(seq_n) ~= 1 || mod(seq_n, 1) ~= 0 || seq_n <= 0
    error('Input arg "seq_n" should be a positive integer!');
end


%% Running mean

running_mean_upd = running_mean + (new_data - running_mean)/seq_n;


return