% AUTHOR: Damir Akchurin
% DATE: 11/14/2022

clear variables
close all
clc

%% INPUT:
% Define the size of the system:
SystemSize = [4, 4];
NumRows = SystemSize(1);
NumCols = SystemSize(2);

% Define the mean of stiffness for each spring:
SpringStiffnessMean = repmat(50, [NumRows, NumCols]);

% Define the standard deviation of stiffness for each spring:
SpringStiffnessSTD = repmat(2.5, [NumRows, NumCols]);

% Define the distribution of stiffness for each spring:
SpringStiffnessDists = repmat("Normal", [NumRows, NumCols]);

% Define the mean of peak strength for each spring:
SpringPeakStrengthMean = repmat(100, [NumRows, NumCols]);

% Define the standard deviation of peak strength for each spring:
SpringPeakStrengthSTD = repmat(10, [NumRows, NumCols]);

% Define the distribution of peak strength for each spring:
SpringPeakStrengthDists = repmat("Normal", [NumRows, NumCols]);

% Define the load applied to the system:
QSystem = 100;

%% ANALYSIS:
SpringStiffnessSamples = zeros(NumRows, NumCols);
SpringPeakStrengthSamples = zeros(NumRows, NumCols);
for i = 1:NumRows
    for j = 1:NumCols
        % Generate stiffness samples for each spring:
        SpringStiffnessSamples(i, j) = random(SpringStiffnessDists(i, j), SpringStiffnessMean(i, j), SpringStiffnessSTD(i, j));

        % Generate peak strength samples for each spring:
        SpringPeakStrengthSamples(i, j) = random(SpringPeakStrengthDists(i, j), SpringPeakStrengthMean(i, j), SpringPeakStrengthSTD(i, j));
    end
end

% Generate stiffness matrices for each spring:
SpringStiffnessMatrices = cell(NumRows, NumCols);
for i = 1:NumRows
    for j = 1:NumCols
        SpringStiffnessMatrices{i, j} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
    end
end

% Generate stiffness matrix for each row of the system:
RowStiffnessMatrices = cell(NumRows, 1);
for i = 1:NumRows
    RowStiffnessMatrices{i} = zeros(2, 2);
    for j = 1:NumCols
        RowStiffnessMatrices{i} = RowStiffnessMatrices{i} + SpringStiffnessMatrices{i, j};
    end
end

% Assemble global stiffness matrix for the whole system:
GlobalStiffnessMatrix = zeros(NumRows + 1);
for i = 1:NumRows
    GlobalStiffnessMatrix(i:i + 1, i:i + 1) = GlobalStiffnessMatrix(i:i + 1, i:i + 1) + RowStiffnessMatrices{i};
end

% Assemble global force vector for the whole system:
GlobalForceVector = zeros(NumRows + 1, 1);
GlobalForceVector(end) = QSystem;

% Apply boundary conditions using the elimination approach:
% 1. Since the top of the system is always fixed:
GlobalForceVector(1) = [];
GlobalStiffnessMatrix(1, :) = [];
GlobalStiffnessMatrix(:, 1) = [];

% Solve for the global displacement vector for the whole system:
GlobalDisplacementVector = [0; GlobalStiffnessMatrix\GlobalForceVector];

% Find the displacement vector for each row of the system:
RowDisplacementVectors = cell(NumRows, 1);
for i = 1:NumRows
    RowDisplacementVectors{i} = GlobalDisplacementVector(i:i + 1);
end

% Find the axial force within each spring:
SpringAxialForces = zeros(NumRows, NumCols);
for i = 1:NumRows
    for j = 1:NumCols
        AxialForceVector = SpringStiffnessMatrices{i, j}*RowDisplacementVectors{i};
        SpringAxialForces(i, j) = AxialForceVector(2);
    end
end
