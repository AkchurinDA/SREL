% AUTHOR: Damir Akchurin
% DATE: 11/25/2022

function [SystemDisplacement, SystemStrength] = AnalyzeSSUsingDCMWithNR(SpringConnectivity, NDOFs,...
    InitialStiffnesses, TPStiffnessCoefficients,...
    TPStrengths, TPDisplacements,...
    DisplacementIncrement,...
    MaxNIncrements, MaxNIterations, Tolerance)

% Compute the number of elements:
NSprings = size(SpringConnectivity, 1);

% Set the reference load vector:
ReferenceLoad = zeros(NDOFs, 1);
ReferenceLoad(1) = -1;
ReferenceLoad(end) = 1;

% Set the initial load factor:
LoadFactorInitial = 0;

% Set the initial total applied load:
TotalAppliedLoadInitial = zeros(NDOFs, 1);

% Set the initial global displacement:
GlobalDisplacementInitial = zeros(NDOFs, 1);

% Set the initial residual force:
ResidualForceInitial = zeros(NDOFs, 1);

% Preallocate:
CurrentElementElongation = cell(MaxNIncrements);
GlobalDisplacementIncrementReference = cell(MaxNIncrements);
GlobalDisplacementIncrementResidual = cell(MaxNIncrements);
LoadFactorIncrement = cell(MaxNIncrements);
LoadFactor = cell(MaxNIncrements, 1);
TotalAppliedLoad = cell(MaxNIncrements, 1);
GlobalDisplacementIncrement = cell(MaxNIncrements, 1);
GlobalDisplacement = cell(MaxNIncrements, 1);
GlobalInternalForce = cell(MaxNIncrements, 1);
ResidualForce = cell(MaxNIncrements, 1);
TotalAppliedLoadRecord = cell(MaxNIncrements, 1);
GlobalDisplacementRecord = cell(MaxNIncrements, 1);

% Start incremental steps:
for Increment = 1:MaxNIncrements
    % Start iterative steps:
    for Iteration = 1:MaxNIterations
        % Display the current status:
        disp(strcat("Increment = ", num2str(Increment), ", Iteration = ", num2str(Iteration)))

        if Iteration == 1
            % Construct the element stiffness matrices:
            CurrentElementStiffness = zeros(NSprings, 1);
            ElementStiffness = cell(NSprings, 1);
            for i = 1:NSprings
                % Compute the current element stiffness:
                if Increment == 1
                    CurrentElementStiffness(i) = InitialStiffnesses(i);
                else
                    CurrentElementStiffness(i) = GetCurrentStiffness(CurrentElementElongation{Increment - 1}{end}(i), TPDisplacements{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));
                end

                % Construct the element stiffness matrices:
                ElementStiffness{i} = CurrentElementStiffness(i)*[1, -1; -1, 1];
            end

            % Construct the global tangent stiffness matrix:
            GlobalStiffness = zeros(NDOFs, NDOFs);
            for i = 1:NSprings
                NI = SpringConnectivity(i, 2);
                NJ = SpringConnectivity(i, 3);

                GlobalStiffness(NI, NI) = GlobalStiffness(NI, NI) + ElementStiffness{i}(1, 1);
                GlobalStiffness(NI, NJ) = GlobalStiffness(NI, NJ) + ElementStiffness{i}(1, 2);
                GlobalStiffness(NJ, NI) = GlobalStiffness(NJ, NI) + ElementStiffness{i}(2, 1);
                GlobalStiffness(NJ, NJ) = GlobalStiffness(NJ, NJ) + ElementStiffness{i}(2, 2);
            end

            % Check for singularity:
            if det(GlobalStiffness(2:end, 2:end)) == 0
                break
            end

            % Compute the global displacement increment due to the reference load:
            GlobalDisplacementIncrementReference{Increment}{Iteration} = [0; GlobalStiffness(2:end, 2:end)\ReferenceLoad(2:end)];

            % Compute the global displacement increment due to the reference load:
            GlobalDisplacementIncrementResidual{Increment}{Iteration} = [0; GlobalStiffness(2:end, 2:end)\ResidualForceInitial(2:end)];

            % Compute the load factor increment parameter:
            LoadFactorIncrement{Increment}{Iteration} = DisplacementIncrement/GlobalDisplacementIncrementReference{Increment}{Iteration}(end);
        else
            % Construct the element stiffness matrices:
            CurrentElementStiffness = zeros(NSprings, 1);
            ElementStiffness = cell(NSprings, 1);
            for i = 1:NSprings
                % Compute the current element stiffness:
                CurrentElementStiffness(i) = GetCurrentStiffness(CurrentElementElongation{Increment}{Iteration - 1}(i), TPDisplacements{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));

                % Construct the element stiffness matrices:
                ElementStiffness{i} = CurrentElementStiffness(i)*[1, -1; -1, 1];
            end

            % Construct the global tangent stiffness matrix:
            GlobalStiffness = zeros(NDOFs, NDOFs);
            for i = 1:NSprings
                NI = SpringConnectivity(i, 2);
                NJ = SpringConnectivity(i, 3);

                GlobalStiffness(NI, NI) = GlobalStiffness(NI, NI) + ElementStiffness{i}(1, 1);
                GlobalStiffness(NI, NJ) = GlobalStiffness(NI, NJ) + ElementStiffness{i}(1, 2);
                GlobalStiffness(NJ, NI) = GlobalStiffness(NJ, NI) + ElementStiffness{i}(2, 1);
                GlobalStiffness(NJ, NJ) = GlobalStiffness(NJ, NJ) + ElementStiffness{i}(2, 2);
            end

            % Check for singularity:
            if det(GlobalStiffness(2:end, 2:end)) == 0
                break
            end

            % Compute the global displacement increment due to the reference load:
            GlobalDisplacementIncrementReference{Increment}{Iteration} = [0; GlobalStiffness(2:end, 2:end)\ReferenceLoad(2:end)];

            % Compute the global displacement increment due to the reference load:
            GlobalDisplacementIncrementResidual{Increment}{Iteration} = [0; GlobalStiffness(2:end, 2:end)\ResidualForce{Increment}{Iteration - 1}(2:end)];

            % Compute the load factor increment parameter:
            LoadFactorIncrement{Increment}{Iteration} = -GlobalDisplacementIncrementResidual{Increment}{Iteration}(end)/GlobalDisplacementIncrementReference{Increment}{Iteration}(end);
        end

        % Update the load factor parameter:
        if Increment == 1 && Iteration == 1
            LoadFactor{Increment}{Iteration} = LoadFactorInitial + LoadFactorIncrement{Increment}{Iteration};
        elseif Increment == 1 && Iteration >= 2
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment}{Iteration - 1} + LoadFactorIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration == 1
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment - 1}{end} + LoadFactorIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration >= 2
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment}{Iteration - 1} + LoadFactorIncrement{Increment}{Iteration};
        end

        % Update the total applied load:
        if Increment == 1 && Iteration == 1
            TotalAppliedLoad{Increment}{Iteration} = TotalAppliedLoadInitial + LoadFactorIncrement{Increment}{Iteration}*ReferenceLoad;
        elseif Increment == 1 && Iteration >= LoadFactorIncrement
            TotalAppliedLoad{Increment}{Iteration} = TotalAppliedLoad{Increment}{Iteration - 1} + LoadFactorIncrement{Increment}{Iteration}*ReferenceLoad;
        elseif Increment >= 2 && Iteration == 1
            TotalAppliedLoad{Increment}{Iteration} = TotalAppliedLoad{Increment - 1}{end} + LoadFactorIncrement{Increment}{Iteration}*ReferenceLoad;
        elseif Increment >= 2 && Iteration >= 2
            TotalAppliedLoad{Increment}{Iteration} = TotalAppliedLoad{Increment}{Iteration - 1} + LoadFactorIncrement{Increment}{Iteration}*ReferenceLoad;
        end

        % Compute the global displacement increment:
        GlobalDisplacementIncrement{Increment}{Iteration} = LoadFactorIncrement{Increment}{Iteration}*GlobalDisplacementIncrementReference{Increment}{Iteration} + GlobalDisplacementIncrementResidual{Increment}{Iteration};

        % Update the total displacement:
        if Increment == 1 && Iteration == 1
            GlobalDisplacement{Increment}{Iteration} = GlobalDisplacementInitial + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment == 1 && Iteration >= 2
            GlobalDisplacement{Increment}{Iteration} = GlobalDisplacement{Increment}{Iteration - 1} + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration == 1
            GlobalDisplacement{Increment}{Iteration} = GlobalDisplacement{Increment - 1}{end} + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration >= 2
            GlobalDisplacement{Increment}{Iteration} = GlobalDisplacement{Increment}{Iteration - 1} + GlobalDisplacementIncrement{Increment}{Iteration};
        end

        % Compute the current element elongation:
        CurrentElementElongation{Increment}{Iteration} = zeros(NSprings, 1);
        for i = 1:NSprings
            NI = SpringConnectivity(i, 2);
            NJ = SpringConnectivity(i, 3);

            CurrentElementElongation{Increment}{Iteration}(i) = GlobalDisplacement{Increment}{Iteration}(NJ) - GlobalDisplacement{Increment}{Iteration}(NI);
        end

        % Compute global internal force:
        GlobalInternalForce{Increment}{Iteration} = zeros(NDOFs, 1);
        for i = 1:NSprings
            NI = SpringConnectivity(i, 2);
            NJ = SpringConnectivity(i, 3);

            ElementInternalForce = GetCurrentInternalForce(CurrentElementElongation{Increment}{Iteration}(i), TPDisplacements{i}, TPStrengths{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));

            GlobalInternalForce{Increment}{Iteration}(NI) = GlobalInternalForce{Increment}{Iteration}(NI) - ElementInternalForce;
            GlobalInternalForce{Increment}{Iteration}(NJ) = GlobalInternalForce{Increment}{Iteration}(NJ) + ElementInternalForce;
        end

        % Compute the unbalanced force:
        ResidualForce{Increment}{Iteration} = TotalAppliedLoad{Increment}{Iteration} - GlobalInternalForce{Increment}{Iteration};

        % Compute the norm of the unbalanced force:
        NormResidualForce = sqrt(dot(ResidualForce{Increment}{Iteration}, ResidualForce{Increment}{Iteration}));

        % Check for convergance:
        if NormResidualForce < Tolerance
            TotalAppliedLoadRecord{Increment} = TotalAppliedLoad{Increment}{Iteration};
            GlobalDisplacementRecord{Increment} = GlobalDisplacement{Increment}{Iteration};
            break
        end
    end

    % Check for singularity:
    if det(GlobalStiffness(2:end, 2:end)) == 0
        disp("Terminated due to singularity of the global stiffness matrix!")
        break
    end
end

% Extract the loads & displacements:
SystemDisplacement = zeros(Increment - 1, 1);
SystemStrength = zeros(Increment - 1, 1);
for i = 1:Increment - 1
    SystemDisplacement(i) = GlobalDisplacementRecord{i}(end);
    SystemStrength(i) = TotalAppliedLoadRecord{i}(end);
end
end