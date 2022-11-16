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
SpringPeakStrengthMean = repmat(30, [NumRows, NumCols]);

% Define the standard deviation of peak strength for each spring:
SpringPeakStrengthSTD = repmat(2.5, [NumRows, NumCols]);

% Define the distribution of peak strength for each spring:
SpringPeakStrengthDists = repmat("Normal", [NumRows, NumCols]);

% Define the displacement applied to the system:
DisplacementSystem = 2;

% Define the number of increments:
NIncrements = 1000;

% Define the number of simulations:
NSimulations = 1E5;

%% ANALYSIS:
% Preallocate variables before looping through the simulations:
SystemStrengthRecord          = zeros(NSimulations, 1);
SystemDisplacementRecord      = zeros(NSimulations, 1);
SystemLoadRecord              = zeros(NSimulations, 1);

% Loop through the simulations:
for s = 1:NSimulations
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

    % Define the displacement increments:
    Displacement = linspace(0, DisplacementSystem, NIncrements);

    % Preallocate variables before looping through the displacement increments:
    SpringStiffnessMatrices         = cell(NumRows, NumCols, NIncrements);
    RowStiffnessMatrices            = cell(NumRows, 1, NIncrements);
    GlobalStiffnessMatrix           = cell(1, NIncrements);
    GlobalForceVector               = cell(1, NIncrements);
    GlobalDisplacementVector        = cell(1, NIncrements);
    RowDisplacementVectors          = cell(NumRows, 1, NIncrements);
    SpringAxialForces               = cell(1, NIncrements);
    FailedSpringsRecord             = zeros(NumRows, NumCols, NIncrements);

    % Loop through the displacement increments:
    for k = 1:NIncrements
        % Generate stiffness matrices for each spring:
        if k == 1
            for i = 1:NumRows
                for j = 1:NumCols
                    SpringStiffnessMatrices{i, j, k} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
                end
            end
        else
            % If any springs failed on the previous displacement iteration, modify their stiffness matrices:
            if any(FailedSpringsRecord(:, :, k - 1), "all")
                for i = 1:NumRows
                    for j = 1:NumCols
                        if FailedSpringsRecord(i, j, k - 1) == 1
                            SpringStiffnessMatrices{i, j, k} = zeros(2, 2);
                        else
                            SpringStiffnessMatrices{i, j, k} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
                        end
                    end
                end
            else
                for i = 1:NumRows
                    for j = 1:NumCols
                        SpringStiffnessMatrices{i, j, k} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
                    end
                end
            end
        end

        % Generate stiffness matrix for each row of the system:
        for i = 1:NumRows
            RowStiffnessMatrices{i, 1, k} = zeros(2, 2);
            for j = 1:NumCols
                RowStiffnessMatrices{i, 1, k} = RowStiffnessMatrices{i, 1, k} + SpringStiffnessMatrices{i, j, k};
            end
        end

        % Assemble global stiffness matrix for the whole system:
        GlobalStiffnessMatrix{k} = zeros(NumRows + 1);
        for i = 1:NumRows
            GlobalStiffnessMatrix{k}(i:i + 1, i:i + 1) = GlobalStiffnessMatrix{k}(i:i + 1, i:i + 1) + RowStiffnessMatrices{i, 1, k};
        end

        % Assemble global force vector for the whole system:
        if k == 1
            GlobalForceVector{k} = zeros(NumRows + 1, 1);
        else
            % If any springs failed on the previous displacement iteration, modify the global force vector:
            if any(FailedSpringsRecord(:, :, k - 1), "all")
                GlobalForceVector{k} = zeros(NumRows + 1, 1);
                for i = 1:NumRows
                    for j = 1:NumCols
                        % Add spring's peak strength as a force to the top row:
                        GlobalForceVector{k}(i) = GlobalForceVector{k}(i) + ...
                            FailedSpringsRecord(i, j, k - 1)*SpringPeakStrengthSamples(i, j);
                        % Add spring's peak strength as a force to the bottom row:
                        GlobalForceVector{k}(i + 1) = GlobalForceVector{k}(i + 1) - ...
                            FailedSpringsRecord(i, j, k - 1)*SpringPeakStrengthSamples(i, j);
                    end
                end
            else
                GlobalForceVector{k} = zeros(NumRows + 1, 1);
            end
        end

        % Apply boundary conditions using the elimination approach:
        % 1. Since the top of the system is always fixed:
        GlobalForceVector{k}(1) = [];
        GlobalStiffnessMatrix{k}(1, :) = [];
        GlobalStiffnessMatrix{k}(:, 1) = [];

        % 2. Since a known displacement is applied at the bottom of the system:
        D = Displacement(k);
        GlobalForceVector{k}(end) = [];
        GlobalStiffnessMatrix{k}(end, :) = [];
        GlobalForceVector{k} = GlobalForceVector{k} - D*GlobalStiffnessMatrix{k}(:, end);
        GlobalStiffnessMatrix{k}(:, end) = [];

        % Solve for the global displacement vector for the whole system:
        GlobalDisplacementVector{k} = [0; GlobalStiffnessMatrix{k}\GlobalForceVector{k}; D];

        % Find the displacement vector for each row of the system:
        for i = 1:NumRows
            RowDisplacementVectors{i, 1, k} = GlobalDisplacementVector{k}(i:i + 1);
        end

        % Find the axial force within each spring:
        SpringAxialForces{k} = zeros(NumRows, NumCols);
        for i = 1:NumRows
            for j = 1:NumCols
                AxialForceVector = SpringStiffnessMatrices{i, j, k}*RowDisplacementVectors{i, 1, k};
                SpringAxialForces{k}(i, j) = AxialForceVector(2);
            end
        end

        % Check if any springs fail:
        for i = 1:NumRows
            for j = 1:NumCols
                % If a spring fails, record it:
                if SpringAxialForces{k}(i, j) >= SpringPeakStrengthSamples(i, j) || (SpringAxialForces{k}(i, j) == 0 && k ~= 1)
                    FailedSpringsRecord(i, j, k) = 1;
                end
            end
        end

        % Check if the systems does fail before reaching the displacement applied to the system:
        if any(sum(FailedSpringsRecord(:, :, k), 2) == NumCols)
            % If a row fails, record it:
            FailedRowIndex = find(sum(FailedSpringsRecord(:, :, k), 2) == NumCols);

            % Find the load applied to the system at failure:
            % 1. Reconstruct the global stiffness matrix:
            GlobalStiffnessMatrixFull = zeros(NumRows + 1);
            for i = 1:NumRows
                GlobalStiffnessMatrixFull(i:i + 1, i:i + 1) = GlobalStiffnessMatrixFull(i:i + 1, i:i + 1) + RowStiffnessMatrices{i, 1, k};
            end

            % 2. Calculate the force applied to the last row of the system:
            GlobalForceVectorFull = GlobalStiffnessMatrixFull*GlobalDisplacementVector{k};
            SystemLoadFailure = GlobalForceVectorFull(end) + sum(FailedSpringsRecord(end, :, k).*SpringPeakStrengthSamples(end, :));

            % Find the displacement of the system at failure:
            SystemDisplacementFailure = D;

            % Find the strength of the system at failure:
            SystemStrengthFailure = min(sum(SpringPeakStrengthSamples(FailedRowIndex, :), 2));
            break
        end

        % Check if the systems does not fail upon reaching the displacement applied to the system:
        if k == NIncrements && any(sum(FailedSpringsRecord(:, :, k), 2) ~= NumCols)
            % Find the load applied to the system:
            % 1. Reconstruct the global stiffness matrix:
            GlobalStiffnessMatrixFull = zeros(NumRows + 1);
            for i = 1:NumRows
                GlobalStiffnessMatrixFull(i:i + 1, i:i + 1) = GlobalStiffnessMatrixFull(i:i + 1, i:i + 1) + RowStiffnessMatrices{i, 1, k};
            end

            % 2. Calculate the force applied to the last row of the system:
            GlobalForceVectorFull = GlobalStiffnessMatrixFull*GlobalDisplacementVector{k};
            SystemLoadFailure = GlobalForceVectorFull(end) + sum(FailedSpringsRecord(end, :, k).*SpringPeakStrengthSamples(end, :));

            % Find the strength of the system:
            % 1. Keep loading the system until it fails to find the strength of the system:
            DisplacementIncrement = DisplacementSystem/NIncrements;
            q = k + 1;
            while any(sum(FailedSpringsRecord(:, :, q - 1), 2) ~= NumCols)
                % Generate stiffness matrices for each spring:
                if any(FailedSpringsRecord(:, :, q - 1), "all")
                    for i = 1:NumRows
                        for j = 1:NumCols
                            if FailedSpringsRecord(i, j, q - 1) == 1
                                SpringStiffnessMatrices{i, j, q} = zeros(2, 2);
                            else
                                SpringStiffnessMatrices{i, j, q} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
                            end
                        end
                    end
                else
                    for i = 1:NumRows
                        for j = 1:NumCols
                            SpringStiffnessMatrices{i, j, q} = SpringStiffnessSamples(i, j)*[1, -1; -1, 1];
                        end
                    end
                end

                % Generate stiffness matrix for each row of the system:
                for i = 1:NumRows
                    RowStiffnessMatrices{i, 1, q} = zeros(2, 2);
                    for j = 1:NumCols
                        RowStiffnessMatrices{i, 1, q} = RowStiffnessMatrices{i, 1, q} + SpringStiffnessMatrices{i, j, q};
                    end
                end

                % Assemble global stiffness matrix for the whole system:
                GlobalStiffnessMatrix{q} = zeros(NumRows + 1);
                for i = 1:NumRows
                    GlobalStiffnessMatrix{q}(i:i + 1, i:i + 1) = GlobalStiffnessMatrix{q}(i:i + 1, i:i + 1) + RowStiffnessMatrices{i, 1, q};
                end

                % Assemble global force vector for the whole system:
                if any(FailedSpringsRecord(:, :, q - 1), "all")
                    GlobalForceVector{q} = zeros(NumRows + 1, 1);
                    for i = 1:NumRows
                        for j = 1:NumCols
                            % Add spring's peak strength as a force to the top row:
                            GlobalForceVector{q}(i) = GlobalForceVector{q}(i) + ...
                                FailedSpringsRecord(i, j, q - 1)*SpringPeakStrengthSamples(i, j);
                            % Add spring's peak strength as a force to the bottom row:
                            GlobalForceVector{q}(i + 1) = GlobalForceVector{q}(i + 1) - ...
                                FailedSpringsRecord(i, j, q - 1)*SpringPeakStrengthSamples(i, j);
                        end
                    end
                else
                    GlobalForceVector{q} = zeros(NumRows + 1, 1);
                end

                % Apply boundary conditions using the elimination approach:
                % 1. Since the top of the system is always fixed:
                GlobalForceVector{q}(1) = [];
                GlobalStiffnessMatrix{q}(1, :) = [];
                GlobalStiffnessMatrix{q}(:, 1) = [];

                % 2. Since a known displacement is applied at the bottom of the system:
                D = DisplacementSystem + DisplacementIncrement*(q - k);
                GlobalForceVector{q}(end) = [];
                GlobalStiffnessMatrix{q}(end, :) = [];
                GlobalForceVector{q} = GlobalForceVector{q} - D*GlobalStiffnessMatrix{q}(:, end);
                GlobalStiffnessMatrix{q}(:, end) = [];

                % Solve for the global displacement vector for the whole system:
                GlobalDisplacementVector{q} = [0; GlobalStiffnessMatrix{q}\GlobalForceVector{q}; D];

                % Find the displacement vector for each row of the system:
                for i = 1:NumRows
                    RowDisplacementVectors{i, 1, q} = GlobalDisplacementVector{q}(i:i + 1);
                end

                % Find the axial force within each spring:
                SpringAxialForces{q} = zeros(NumRows, NumCols);
                for i = 1:NumRows
                    for j = 1:NumCols
                        AxialForceVector = SpringStiffnessMatrices{i, j, q}*RowDisplacementVectors{i, 1, q};
                        SpringAxialForces{q}(i, j) = AxialForceVector(2);
                    end
                end

                % Check if any springs fail:
                FailedSpringsRecord(:, :, q) = zeros(NumRows, NumCols);
                for i = 1:NumRows
                    for j = 1:NumCols
                        % If a spring fails, record it:
                        if SpringAxialForces{q}(i, j) >= SpringPeakStrengthSamples(i, j) || SpringAxialForces{q}(i, j) == 0
                            FailedSpringsRecord(i, j, q) = 1;
                        end
                    end
                end

                % Check if the systems does fail before reaching the displacement applied to the system:
                if any(sum(FailedSpringsRecord(:, :, q), 2) == NumCols)
                    % If a row fails, record it:
                    FailedRowIndex = find(sum(FailedSpringsRecord(:, :, q), 2) == NumCols);

                    % Find the strength of the system at failure:
                    SystemStrengthFailure = min(sum(SpringPeakStrengthSamples(FailedRowIndex, :), 2));
                    break
                else
                    q = q + 1;
                end
            end

            % Find the displacement of the system:
            SystemDisplacementFailure = D;
            break
        end
    end

    % Store values for each simulation:
    SystemStrengthRecord(s) = SystemStrengthFailure;
    SystemDisplacementRecord(s) = SystemDisplacementFailure;
    SystemLoadRecord(s) = SystemLoadFailure;
end

% Store all values:
save(strcat("HSS_FEM_DC_", num2spr(NSimulations)), "SystemStrengthRecord", "SystemDisplacementFailure", "SystemLoadRecord")

