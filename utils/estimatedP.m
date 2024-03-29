function estP = estimatedP(realV, randomV, tailed)
%% Probability estimation given a set of null model values
%
% USAGE: estP = estimatedP(realV, randomVs, tailed=2)
%
% Given a real value of some statistic / feature and a generated set of 
% null-model values, the function estimates the probability of 
% the real value. 
% 
% Works on multiple results as well - if there are N tests to perform,
% realV could be a vector and randomV a matrix with N columns.
%
% Inputs:
% realV         -  Vector of length N, numeric. Contains the values we
%               compare against the random values generated from some 
%               null-model
% randomV       - Matrix with N columns, each column j containing the random
%               values generated by a null-model corresponding to realV(j)
% tailed        - One vs. two-tailed test, plus the direction of our hypothesis. 
%               Numeric, one of [-1, 1, 2]. -1 is a left-tail test meaning 
%               that the null hypothesis is that realV >= randV and that we 
%               want to reject that. +1 means the opposite, right-tailed 
%               expectation; +2 means two-tailed test with appropriate 
%               correction of the p-value. Defaults to 2.
%
% Output:
% estP          - Vector of length N containing the estimated p values,
%               numeric .If arg tailed==2, there is also a second column of
%               values specifying the direction of the difference.
%               Direction is coded by -1 and +1, corresponding to
%               realV<randomV and realV>randV, respectively
%
%
% Notes:
% (1) We know about but do not follow Phipson and Smith
% (http://www.statsci.org/smyth/pubs/permp.pdf) and do let p=0 values as
% the result
% (2) We treat equality (realV=randV), but know that it practically never 
% happens with float numerics and in our use cases
%


%% Input checks

% numbe rof inputs, assign default value to arg tailed if necessary
if ~ismember(nargin, [2 3])
    error('Need input args "realV" and "randomV", also optional arg "tailed"!');
end
if nargin == 2
    tailed = 2;
end
% realV and randomV should a vector and a matrix
if ~isvector(realV)
    error('Input arg "realV" should be a vector!');
end
if ~ismatrix(randomV)
    error('Input arg "randomV" should be a matrix!');
end
% tailed is one of [-1, 1, 2]
if ~ismember(tailed, [-1, 1 ,2])
    error('Input arg "tailed" should be one of [-1 ,1, 2]!');
end
% provide warning if there are more tests then random values - maybe randV
% should be transposed?
if size(randomV, 1) <= size(randomV, 2)
    warning('There seem to be more tests than random values - perhaps input arg "randV" should be transposed?');
end
% transpose realV if it is a column vector
if size(realV, 1) > size(realV, 2)
    realV = realV';
end
    

%% Estimate

% to do the tests properly with regards to equality (equality should always 
% boost the null hypothesis) we handle each case of arg tailed separately 

switch tailed
    case -1  % what is the probability that the real value is not smaller than the random values?
        estP = sum(randomV<=realV)/size(randomV, 1);
    case 1  % what is the probability that the real value is not larger than the random values?
        estP = sum(randomV>=realV)/size(randomV, 1);
    case 2  % do both, choose the smaller value, do correction
        estProbSmaller = sum(randomV<=realV)/size(randomV, 1);
        estProbLarger = sum(randomV>=realV)/size(randomV, 1);
        estP(:, 1) = min(estProbSmaller, estProbLarger)*2;  % with correction
        % second column tells us about the direction of the result
        estP(:, 2) = ones(length(realV), 1);
        estP(estProbSmaller<estProbLarger, 2) = -1;
end


return













