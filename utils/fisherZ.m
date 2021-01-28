function transformedArray = fisherZ(corrArray)
%% Fisher-Z transformation for correlation values
%
% Placeholder and reminder, for those who - like me - keep forgetting 
% that Fisher-z is simply atanh(). Though the formulation used here is
% slightly faster...
% 
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
if ~isnumeric(corrArray) || any(corrArray(:)>1) || any(corrArray(:)<-1)
    error('Input arg "corrArray" must a numeric array containing values only in the range [-1, 1] (correlation values)!');
end

%% Transform

% The classic Fisher-Z formulation:
transformedArray = (log((1+corrArray)./(1-corrArray)))/2;

% % Simpler solution in matlab:
% transformedArray = atanh(corrArray);


return
