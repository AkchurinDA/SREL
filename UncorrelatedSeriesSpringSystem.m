% Author: Damir Akchurin
% Date: 1/1/22
% Version history:
% V1.0 - First working version of the software (1/2/22)

clc
close all
clear variables

% Define the number of springs:
NSprings = 10;

% Define the properties of the springs:
ReliabilityIndexMin = 2.5;
ReliabilityIndexMax = 4.0;
ReliabilityIndices = ReliabilityIndexMin + (ReliabilityIndexMax - ReliabilityIndexMin)*rand(NSprings, 1);
FailureProbabilies = normcdf(-ReliabilityIndices); % Assuming that the strength distribution is normal for each spring

% Calculate the system failure probability:
SystemFailureProbability = 1 - prod(1 - FailureProbabilies);
SystemReliabilityIndex = -norminv(SystemFailureProbability); % Assuming that the strength distribution is normal for the whole system