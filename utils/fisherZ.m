function transformedArray = fisherZ(corrArray)
%% Fisher-Z transformation for correlation values
%
% USAGE: transformedArray = fisherZ(corrArray)
%
% Input:
% corrArray         - Numeric array of correlation values. Values must be
%                   in range [-1 1].
%
% Output:
% transformedArray  - Numeric array, Fisher z-transformed correlation
%                   values. Same size as input arg "corrArray".
%


%% Input checks

if nargin ~= 1
    error('Function fisherZ requires input arg "corrArray"!');
end
tooLarge = corrArray>1; tooSmall = corrArray<-1;
if ~isnumeric(corrArray) || any(tooLarge(:)) || any(tooSmall(:))
    error('Input arg "corrArray" must a numeric array containing values only in the range [-1, 1] (correlation values)!');
end

%% Transform

transformedArray = (log((1+corrArray)./(1-corrArray)))/2;



return
