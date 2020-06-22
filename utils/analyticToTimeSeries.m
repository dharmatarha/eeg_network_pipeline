function [realData] = analyticToTimeSeries(angleData, envData)
%% Calculate real-valued signal from the angle and envelope (magnitude) of the analytic signal 
% 
% CONSIDER CALLING BUILT-IN pol2cart DIRECTLY
% 
% USAGE: [realData] = analyticToTimeSeries(angleData, envData)
%
% Calculates real-valued signal from polar representation of the analytic
% signal (angle in radians and magnitude / envelope). 
% Works with any matrix/tensor. Simply calls pol2cart.
%
%
% Input: 
% angleData     - Numerical matrix (real), angle data in radians. Same size
%               as "envData"
% envData       - Numerical matrix (real), magnitude data. Same size as 
%               "angleData"
%
% Output:
% realData      - Numerical matrix (real), same size as "angleData" and 
%               "envData". Real-valued part (original time series) of the 
%               analytic signal.

[realData, ~] = pol2cart(angleData, envData);

return

