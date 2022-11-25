% AUTHOR: Damir Akchurin
% DATE: 11/25/2022

function [SystemDisplacement, SystemStrength] = AnalyzeSSUsingGDCM(SpringConnectivity, NDOFs,...
    InitialStiffnesses, TPStiffnessCoefficients,...
    TPStrengths, TPDisplacements,...
    MaxNIncrements, MaxNIterations, Tolerance)

% Compute the number of elements:
NSprings = size(SpringConnectivity, 1);

% Set the reference load vector:
ReferenceLoad = zeros(NDOFs, 1);
ReferenceLoad(1) = -1;
ReferenceLoad(end) = 1;

% Set the maximum load factor:
MaxLoadFactor = 1000;

% Set the initial load factor:
LoadFactorInitial = 0;

% Set the initial load increment:
LoadIncrementInitial = 1;

% Set the initial total displacement:
TotalDisplacementInitial = zeros(NDOFs, 1);

% Preallocate memory for the global displacement increment due to the reference load:
GlobalDisplacementIncrementReference = cell(MaxNIncrements, 1);

% Preallocate memory for the global displacement increment due to the residual load:
GlobalDisplacementIncrementResidual = cell(MaxNIncrements, 1);

% Preallocate memory for the load increment parameter:
LoadIncrement = cell(MaxNIncrements, 1);

% Preallocate memory for the global displacement increment:
GlobalDisplacementIncrement = cell(MaxNIncrements, 1);

% Preallocate memory for the load factor parameter:
LoadFactor = cell(MaxNIncrements, 1);

% Preallocate memory for the total displacement:
TotalDisplacement = cell(MaxNIncrements, 1);

% Preallocate memory for the total applied load:
TotalAppliedLoad = cell(MaxNIncrements, 1);

% Preallocate memory for the global internal force:
GlobalInternalForce = cell(MaxNIncrements, 1);

% Preallocate memory for the unbalanced force:
ResidualForce = cell(MaxNIncrements, 1);

% Start incremental steps:
for Increment = 1:MaxNIncrements
    % Set the initial residual force to zero since each increment step start from the equilibrium state:
    ResidualForceInitial = zeros(NDOFs, 1);

    % Start iterative steps:
    for Iteration = 1:MaxNIterations
        % Displat current status:
        disp(strcat("Increment = ", num2str(Increment), ", Iteration = ", num2str(Iteration)))

        % First iteration:
        if Iteration == 1
            % Construct the element stiffness matrices:
            CurrentElementStiffness = zeros(NSprings, 1);
            ElementStiffness = cell(NSprings, 1);
            for i = 1:NSprings
                if Increment == 1
                    CurrentElementStiffness(i) = InitialStiffnesses(i);
                else
                    CurrentElementStiffness(i) = GetCurrentStiffness(CurrentElementElongation(i), TPDisplacements{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));
                end
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

            % Compute the load increment parameter:
            if Increment == 1
                GSP = 1;
            else
                GSP = dot(GlobalDisplacementIncrementReference{1}{1}, GlobalDisplacementIncrementReference{1}{1})/dot(GlobalDisplacementIncrementReference{Increment - 1}{1}, GlobalDisplacementIncrementReference{Increment}{1});
            end

            if GSP > 0
                LoadIncrement{Increment}{Iteration} = LoadIncrementInitial*sqrt(abs(GSP));
            else
                LoadIncrement{Increment}{Iteration} = -LoadIncrementInitial*sqrt(abs(GSP));
            end
        else
            % Construct the element stiffness matrices:
            CurrentElementStiffness = zeros(NSprings, 1);
            ElementStiffness = cell(NSprings, 1);
            for i = 1:NSprings
                if Increment == 1
                    CurrentElementStiffness(i) = InitialStiffnesses(i);
                else
                    CurrentElementStiffness(i) = GetCurrentStiffness(CurrentElementElongation(i), TPDisplacements{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));
                end
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

            % Compute the load increment parameter:
            LoadIncrement{Increment}{Iteration} = -dot(GlobalDisplacementIncrementReference{Increment - 1}{1}, GlobalDisplacementIncrementResidual{Increment}{Iteration})/dot(GlobalDisplacementIncrementReference{Increment - 1}{1}, GlobalDisplacementIncrementReference{Increment}{Iteration});
        end

        % Compute the global displacement increment:
        GlobalDisplacementIncrement{Increment}{Iteration} = LoadIncrement{Increment}{Iteration}*GlobalDisplacementIncrementReference{Increment}{Iteration} + GlobalDisplacementIncrementResidual{Increment}{Iteration};

        % Compute the load factor parameter:
        if Increment == 1 && Iteration == 1
            LoadFactor{Increment}{Iteration} = LoadFactorInitial + LoadIncrement{Increment}{Iteration};
        elseif Increment == 1 && Iteration >= 2
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment}{Iteration - 1} + LoadIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration == 1
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment - 1}{end} + LoadIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration >= 2
            LoadFactor{Increment}{Iteration} = LoadFactor{Increment}{Iteration - 1} + LoadIncrement{Increment}{Iteration};
        end

%         if Iteration == 1
%             LoadFactor{Increment}{Iteration} = LoadFactorInitial + LoadIncrement{Increment}{Iteration};
%         else
%             LoadFactor{Increment}{Iteration} = LoadFactor{Increment}{Iteration - 1} + LoadIncrement{Increment}{Iteration};
%         end

        % Update the total displacement:
        if Increment == 1 && Iteration == 1
            TotalDisplacement{Increment}{Iteration} = TotalDisplacementInitial + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment == 1 && Iteration >= 2
            TotalDisplacement{Increment}{Iteration} = TotalDisplacement{Increment}{Iteration - 1} + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration == 1
            TotalDisplacement{Increment}{Iteration} = TotalDisplacement{Increment - 1}{end} + GlobalDisplacementIncrement{Increment}{Iteration};
        elseif Increment >= 2 && Iteration >= 2
            TotalDisplacement{Increment}{Iteration} = TotalDisplacement{Increment}{Iteration - 1} + GlobalDisplacementIncrement{Increment}{Iteration};
        end

%         if Iteration == 1
%             TotalDisplacement{Increment}{Iteration} = TotalDisplacementInitial + GlobalDisplacementIncrement{Increment}{Iteration};
%         else
%             TotalDisplacement{Increment}{Iteration} = TotalDisplacement{Increment}{Iteration - 1} + GlobalDisplacementIncrement{Increment}{Iteration};
%         end

        % Update the total applied load:
        TotalAppliedLoad{Increment}{Iteration} = LoadFactor{Increment}{Iteration}*ReferenceLoad;

        % Compute global internal force:
        CurrentElementElongation = zeros(NSprings, 1);
        GlobalInternalForce{Increment}{Iteration} = zeros(NDOFs, 1);
        for i = 1:NSprings
            NI = SpringConnectivity(i, 2);
            NJ = SpringConnectivity(i, 3);

            CurrentElementElongation(i) = TotalDisplacement{Increment}{Iteration}(NJ) - TotalDisplacement{Increment}{Iteration}(NI);
            ElementInternalForce = GetCurrentInternalForce(CurrentElementElongation(i), TPDisplacements{i}, TPStrengths{i}, TPStiffnessCoefficients(i, :), InitialStiffnesses(i));

            GlobalInternalForce{Increment}{Iteration}(NI) = GlobalInternalForce{Increment}{Iteration}(NI) - ElementInternalForce;
            GlobalInternalForce{Increment}{Iteration}(NJ) = GlobalInternalForce{Increment}{Iteration}(NJ) + ElementInternalForce;
        end

        disp(CurrentElementElongation)

        % Compute the unbalanced force:
        ResidualForce{Increment}{Iteration} = TotalAppliedLoad{Increment}{Iteration} - GlobalInternalForce{Increment}{Iteration};

        % Compute the norm of the unbalanced force:
        NormResidualForce = sqrt(dot(ResidualForce{Increment}{Iteration}, ResidualForce{Increment}{Iteration}));

        % Check for convergance:
        if NormResidualForce < Tolerance
            break
        end
    end

    % Check for termination:
    if det(GlobalStiffness(2:end, 2:end)) == 0
        break
    end
    
    if LoadFactor{Increment}{Iteration} > MaxLoadFactor
        disp("Terminated!")
        break
    end
end

% Extract the loads & displacements:
SystemDisplacement = zeros(Increment - 1, 1);
SystemStrength = zeros(Increment - 1, 1);
for i = 1:Increment - 1
    SystemDisplacement(i) = TotalDisplacement{i}{1}(end);
    SystemStrength(i) = TotalAppliedLoad{i}{1}(end);
end
end