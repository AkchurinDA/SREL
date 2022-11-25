function [SystemSupportedLoadRecord, SystemDisplacementRecord, SpringInflectionDisplacements, SpringInflectionStrength, NumSprings] = AnalyzeSSSecant(SpringConnectivity, NumDOFs, SpringStiffnessMeans, SpringStiffnessSTDs, SpringStiffnessDistributions, SpringNumInflectionPoints, SpringInflectionStrengthMeans, SpringInflectionStrengthSTDs, SpringInflectionStrengthDistributions, SpringStiffnessCoefficients, DisplacementIncrement, NumIterations)
% Calculate the number of springs in the system:
NumSprings = size(SpringConnectivity, 1);

% Generate samples of stiffness and strengths for each spring:
SpringStiffness = zeros(NumSprings, 1);
SpringInflectionStrength = cell(NumSprings, 1);
for i = 1:NumSprings
    SpringStiffness(i) = random(SpringStiffnessDistributions(i, 2), SpringStiffnessMeans(i, 2), SpringStiffnessSTDs(i, 2));
    for j = 1:SpringNumInflectionPoints(i, 2)
        SpringInflectionStrength{i}(end + 1) = random(SpringInflectionStrengthDistributions(i, j + 1), SpringInflectionStrengthMeans(i, j + 1), SpringInflectionStrengthSTDs(i, j + 1));
    end
end

% Calculate the inflection displacements for each spring:
SpringInflectionDisplacements = cell(NumSprings, 1);
for i = 1:NumSprings
    for j = 1:SpringNumInflectionPoints(i, 2)
        if j == 1
            SpringInflectionDisplacements{i} = SpringInflectionStrength{i}(j)/(SpringStiffnessCoefficients(i, j + 1)*SpringStiffness(i));
        else
            SpringInflectionDisplacements{i}(end + 1) = SpringInflectionDisplacements{i}(end) + (SpringInflectionStrength{i}(j) - SpringInflectionStrength{i}(j - 1))/(SpringStiffnessCoefficients(i, j + 1)*SpringStiffness(i));
        end
    end
end

% Initilize initial state for each spring:
SpringCurrentState = ones(NumSprings, 1);

% Initialize the waitbar:
% WaitbarHandle = waitbar(0, "Running analysis");

% Incrementally apply the displacement to the system:
Iteration = 1; % Initialize the iteration counter
SystemSupportedLoadRecord = zeros(NumIterations, 1);
SystemDisplacementRecord = zeros(NumIterations, 1);
while Iteration <= NumIterations
    if Iteration == 1
        % Calculate the displacement at the current iteration:
        CurrentDisplacement = DisplacementIncrement*Iteration;

        % Generate stiffness matrix for each spring:
        SpringStiffnessMatrix = cell(NumSprings, 1);
        for i = 1:NumSprings
            SpringStiffnessMatrix{i} = SpringStiffness(i)*[1, -1; -1, 1];
        end

        % Assemble the global stiffness matrix:
        GlobalStiffnessMatrix = AssembleGlobalStiffnessMatrix(SpringConnectivity, SpringStiffnessMatrix, NumSprings, NumDOFs);

        % Assemble the global force vector:
        GlobalForceVector = zeros(NumDOFs, 1);

        % Apply the boundary condition to the global stiffness matrix and global force vector:
        GlobalStiffnessMatrixBC = GlobalStiffnessMatrix;
        GlobalForceVectorBC = GlobalForceVector;
        GlobalForceVectorBC(1) = [];
        GlobalStiffnessMatrixBC(1, :) = [];
        GlobalStiffnessMatrixBC(:, 1) = [];
        GlobalForceVectorBC(end) = [];
        GlobalForceVectorBC = GlobalForceVectorBC - CurrentDisplacement*GlobalStiffnessMatrixBC(1:end - 1, end);
        GlobalStiffnessMatrixBC(end, :) = [];
        GlobalStiffnessMatrixBC(:, end) = [];

        % Calculate the global displacement vector:
        GlobalDisplacementVector = [0; GlobalStiffnessMatrixBC\GlobalForceVectorBC; CurrentDisplacement];

        % Calculate the local displacement vector for each spring:
        SpringLocalDisplacementVector = cell(NumSprings, 1);
        for i = 1:NumSprings
            NI = SpringConnectivity(i, 2);
            NJ = SpringConnectivity(i, 3);

            SpringLocalDisplacementVector{i} = [GlobalDisplacementVector(NI); GlobalDisplacementVector(NJ)];
        end

        % Calculate the elongation for each spring:
        SpringElongation = zeros(NumSprings, 1);
        for i = 1:NumSprings
            SpringElongation(i) = SpringLocalDisplacementVector{i}(2) - SpringLocalDisplacementVector{i}(1);
        end

        % Calculate the axial force for each spring:
        SpringAxialForce = zeros(NumSprings, 1);
        for i = 1:NumSprings
            SpringAxialForceVector = SpringStiffnessMatrix{i}*SpringLocalDisplacementVector{i};
            SpringAxialForce(i) = SpringAxialForceVector(2);
        end
    else
        % Calculate the displacement at the current iteration:
        CurrentDisplacement = DisplacementIncrement*Iteration;

        % Calculate current stiffness for each spring:
        CurrentSpringStiffness = zeros(NumSprings, 1);
        for i = 1:NumSprings
            CurrentSpringStiffness(i) = GetCurrentSpringStiffness(SpringElongation(i), SpringInflectionDisplacements{i}, SpringStiffnessCoefficients(i, 2:end), SpringStiffness(i));
        end

        % Calculate current axial force for each spring:
        CurrentSpringAxialForce = zeros(NumSprings, 1);
        for i = 1:NumSprings
            CurrentSpringAxialForce(i) = GetCurrentSpringAxialForce(SpringElongation(i), CurrentSpringStiffness(i), SpringInflectionDisplacements{i}, SpringInflectionStrength{i});
        end

        % Generate stiffness matrix for each spring:
        CurrentSpringStiffness = zeros(NumSprings, 1);
        SpringStiffnessMatrix = cell(NumSprings, 1);
        for i = 1:NumSprings
            CurrentSpringStiffness(i) = CurrentSpringAxialForce(i)/SpringElongation(i);
            SpringStiffnessMatrix{i} = CurrentSpringStiffness(i)*[1, -1; -1, 1];
        end

        % Assemble the global stiffness matrix:
        GlobalStiffnessMatrix = AssembleGlobalStiffnessMatrix(SpringConnectivity, SpringStiffnessMatrix, NumSprings, NumDOFs);

        % Assemble the global force vector:
        GlobalForceVector = zeros(NumDOFs, 1);

        % Apply the boundary condition to the global stiffness matrix and global force vector:
        GlobalStiffnessMatrixBC = GlobalStiffnessMatrix;
        GlobalForceVectorBC = GlobalForceVector;
        GlobalForceVectorBC(1) = [];
        GlobalStiffnessMatrixBC(1, :) = [];
        GlobalStiffnessMatrixBC(:, 1) = [];
        GlobalForceVectorBC(end) = [];
        GlobalForceVectorBC = GlobalForceVectorBC - CurrentDisplacement*GlobalStiffnessMatrixBC(1:end - 1, end);
        GlobalStiffnessMatrixBC(end, :) = [];
        GlobalStiffnessMatrixBC(:, end) = [];

        % Calculate the global displacement vector:
        GlobalDisplacementVector = [0; GlobalStiffnessMatrixBC\GlobalForceVectorBC; CurrentDisplacement];

        % Calculate the local displacement vector for each spring:
        SpringLocalDisplacementVector = cell(NumSprings, 1);
        for i = 1:NumSprings
            NI = SpringConnectivity(i, 2);
            NJ = SpringConnectivity(i, 3);

            SpringLocalDisplacementVector{i} = [GlobalDisplacementVector(NI); GlobalDisplacementVector(NJ)];
        end

        % Calculate the elongation for each spring:
        SpringElongation = zeros(NumSprings, 1);
        for i = 1:NumSprings
            SpringElongation(i) = SpringLocalDisplacementVector{i}(2) - SpringLocalDisplacementVector{i}(1);
        end

        % Calculate the axial force for each spring:
        SpringAxialForce = zeros(NumSprings, 1);
        for i = 1:NumSprings
            SpringAxialForceVector = SpringStiffnessMatrix{i}*SpringLocalDisplacementVector{i};
            SpringAxialForce(i) = SpringAxialForceVector(2);
        end

        % Check if there are any unbalanced forces:
        UnbalancedForceSwitch = zeros(NumSprings, 1);
        for i = 1:NumSprings
            if SpringCurrentState(i) < SpringNumInflectionPoints(i, 2) + 1
                if any(SpringElongation(i) > SpringInflectionDisplacements{i}(SpringCurrentState(i)))
                    UnbalancedForceSwitch(i) = 1;
                    SpringCurrentState(i) = min(SpringCurrentState(i) + 1, SpringNumInflectionPoints(i, 2) + 1);
                end
            end
        end

        % If there are any unbalanced forces within the system, try to reequilibrate them:
        if any(UnbalancedForceSwitch)
            while any(UnbalancedForceSwitch)
                % Calculate current stiffness for each spring:
                CurrentSpringStiffness = zeros(NumSprings, 1);
                for i = 1:NumSprings
                    CurrentSpringStiffness(i) = GetCurrentSpringStiffness(SpringElongation(i), SpringInflectionDisplacements{i}, SpringStiffnessCoefficients(i, 2:end), SpringStiffness(i));
                end

                % Calculate current axial force for each spring:
                CurrentSpringAxialForce = zeros(NumSprings, 1);
                for i = 1:NumSprings
                    CurrentSpringAxialForce(i) = GetCurrentSpringAxialForce(SpringElongation(i), CurrentSpringStiffness(i), SpringInflectionDisplacements{i}, SpringInflectionStrength{i});
                end

                % Generate stiffness matrix for each spring:
                CurrentSpringStiffness = zeros(NumSprings, 1);
                SpringStiffnessMatrix = cell(NumSprings, 1);
                for i = 1:NumSprings
                    CurrentSpringStiffness(i) = CurrentSpringAxialForce(i)/SpringElongation(i);
                    SpringStiffnessMatrix{i} = CurrentSpringStiffness(i)*[1, -1; -1, 1];
                end

                % Assemble the global stiffness matrix:
                GlobalStiffnessMatrix = AssembleGlobalStiffnessMatrix(SpringConnectivity, SpringStiffnessMatrix, NumSprings, NumDOFs);

                % Assemble the global force vector:
                GlobalForceVector = zeros(NumDOFs, 1);

                % Apply the boundary condition to the global stiffness matrix and global force vector:
                GlobalStiffnessMatrixBC = GlobalStiffnessMatrix;
                GlobalForceVectorBC = GlobalForceVector;
                GlobalForceVectorBC(1) = [];
                GlobalStiffnessMatrixBC(1, :) = [];
                GlobalStiffnessMatrixBC(:, 1) = [];
                GlobalForceVectorBC(end) = [];
                GlobalForceVectorBC = GlobalForceVectorBC - CurrentDisplacement*GlobalStiffnessMatrixBC(1:end - 1, end);
                GlobalStiffnessMatrixBC(end, :) = [];
                GlobalStiffnessMatrixBC(:, end) = [];

                % Calculate the global displacement vector:
                GlobalDisplacementVector = [0; GlobalStiffnessMatrixBC\GlobalForceVectorBC; CurrentDisplacement];

                % Calculate the equivalent force appied to the system:
                SystemEquivalentLoad = GlobalStiffnessMatrix(end, :)*GlobalDisplacementVector;

                % Calculate the local displacement vector for each spring:
                SpringLocalDisplacementVector = cell(NumSprings, 1);
                for i = 1:NumSprings
                    NI = SpringConnectivity(i, 2);
                    NJ = SpringConnectivity(i, 3);

                    SpringLocalDisplacementVector{i} = [GlobalDisplacementVector(NI); GlobalDisplacementVector(NJ)];
                end

                % Calculate the elongation for each spring:
                SpringElongation = zeros(NumSprings, 1);
                for i = 1:NumSprings
                    SpringElongation(i) = SpringLocalDisplacementVector{i}(2) - SpringLocalDisplacementVector{i}(1);
                end

                % Calculate the axial force for each spring:
                SpringAxialForce = zeros(NumSprings, 1);
                for i = 1:NumSprings
                    SpringAxialForceVector = SpringStiffnessMatrix{i}*SpringLocalDisplacementVector{i};
                    SpringAxialForce(i) = SpringAxialForceVector(2);
                end

                % Check if there are any unbalanced forces:
                UnbalancedForceSwitch = zeros(NumSprings, 1);
                for i = 1:NumSprings
                    if SpringCurrentState(i) < SpringNumInflectionPoints(i, 2) + 1
                        if any(SpringElongation(i) > SpringInflectionDisplacements{i}(SpringCurrentState(i)))
                            UnbalancedForceSwitch(i) = 1;
                            SpringCurrentState(i) = min(SpringCurrentState(i) + 1, SpringNumInflectionPoints(i, 2) + 1);
                        end
                    end
                end
            end
        end
    end

    % Calculate the force supported by the system:
    TopSpringsIndex = SpringConnectivity(:, 2) == 1;
    SystemSupportedLoad = sum(SpringAxialForce(TopSpringsIndex));

    % Record the force supported by the system:
    SystemSupportedLoadRecord(Iteration) = SystemSupportedLoad;

    % Record the displacements applied to the system:
    SystemDisplacementRecord(Iteration) = CurrentDisplacement;

    % Increase the iteration counter:
    Iteration = Iteration + 1;
end

% Adjust the values:
SystemSupportedLoadRecord = [0; SystemSupportedLoadRecord];
SystemDisplacementRecord = [0; SystemDisplacementRecord];
SystemSupportedLoadRecord(end) = [];
SystemDisplacementRecord(end) = [];
end