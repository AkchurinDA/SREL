% AUTHOR: Damir Akchurin
% DATE: 11/25/2022

clear variables
close all
clc

%% Input:
% Define the spring connectivity in the system:
% NOTE: The top node of the system should always have DOF #1
% NOTE: The bottom node of the system should always have the largest DOF #
SpringConnectivity = [
    1   1   2;
    2   1   2;
    3   1   2;
    4   1   2;
    5   2   3;
    6   2   3;
    7   2   3;
    8   2   3;
    9   3   4;
    10  3   4;
    11  3   4;
    12  3   4;
    13  4   5;
    14  4   5;
    15  4   5;
    16  4   5];

NSprings = size(SpringConnectivity, 1);

% Define the total number of degrees of freedom in the system:
NDOFs = 5;

% Define the mean of stiffness for each spring:
InitialStiffnessMeans = [(1:NSprings)', repmat(45, [NSprings, 1])];

% Define the standard deviation of stiffness for each spring:
InitialStiffnessSTDs = [(1:NSprings)', 5*ones(NSprings, 1)];

% Define the distribution of stiffness for each spring:
InitialStiffnessDistribution = [(1:NSprings)', repmat("Normal", [NSprings, 1])];

% Define the stiffness coefficients for each spring:
TPStiffnessCoefficients = [(1:NSprings)', ones(NSprings, 1), zeros(NSprings, 1)];

% Define the numxber of inflection points for each spring:
NTurnPoints = [(1:NSprings)', 1*ones(NSprings, 1)];

% Define the means of inflection strengths for each spring:
TPStrengthDistributionsMeans = [(1:NSprings)', 25*ones(NSprings, 1)];

% Define the standard deviations of inflection strengths for each spring:
TPStrengthDistributionsSTDs = [(1:NSprings)', 5*ones(NSprings, 1)];

% Define the distribution of inflection strengths for each spring:
TPStrengthDistributions = [(1:NSprings)', repmat("Normal", [NSprings, 1])];

% Define the displacement increment to apply to the system:
DisplacementIncrement = 0.2;

% Define the number of increments:
MaxNIncrements = 50;

% Define the number of iterations:
MaxNIterations = 500;

% Define the tolerance for convergance:
Tolerance = 10E-10;

% Define the mean of the system load:
SystemLoadMean = 70;

% Define the standard deviation of the system load:
SystemLoadSTD = 5;

% Define the distribution of the system load:
SystemLoadDistribution = "Normal";

% Define the number of simulations to run:
NumSimulations = 500;

%% Analyze:
% Simulate the system:
UltimateSystemStrength = zeros(NumSimulations, 1);
SystemLoad = zeros(NumSimulations, 1); 
for k = 1:NumSimulations
    % Generate samples of initial stiffness:
    InitialStiffnesses = zeros(NSprings, 1);
    for i = 1:NSprings
        InitialStiffnesses(i) = random(InitialStiffnessDistribution(i, 2), InitialStiffnessMeans(i, 2), InitialStiffnessSTDs(i, 2));
    end

    % Sample the turn-point strengths and compute the corresponding turn-point displacements:
    TPStrengths = cell(NSprings, 1);
    TPDisplacements = cell(NSprings, 1);
    for i = 1:NSprings
        for j = 1:NTurnPoints(i, 2)
            TPStrengths{i}(end + 1) = random(TPStrengthDistributions(i, j + 1), TPStrengthDistributionsMeans(i, j + 1), TPStrengthDistributionsSTDs(i, j + 1));

            if j == 1
                TPDisplacements{i} = TPStrengths{i}(j)/(TPStiffnessCoefficients(i, j + 1)*InitialStiffnesses(i));
            else
                TPDisplacements{i}(end + 1) = TPDisplacements{i}(end) + (TPStrengths{i}(j) - TPStrengths{i}(j - 1))/(TPStiffnessCoefficients(i, j + 1)*InitialStiffnesses(i));
            end
        end
    end

    % Analyze the system:
    [SystemDisplacement, SystemStrength] = AnalyzeSSUsingDCMWithMNR(SpringConnectivity, NDOFs, InitialStiffnesses, TPStiffnessCoefficients, TPStrengths, TPDisplacements, DisplacementIncrement, MaxNIncrements, MaxNIterations, Tolerance);

    % Extract the ultimate strength of the system:
    UltimateSystemStrength(k) = max(SystemStrength);

    % Generate samples of the system load:
    SystemLoad(k) = random(SystemLoadDistribution, SystemLoadMean, SystemLoadSTD);
end

% Compute the limit state function:
LimitStateFunction = UltimateSystemStrength - SystemLoad;

% Compute the probability of failure:
FailureProbability = numel(LimitStateFunction(LimitStateFunction <= 0))/numel(LimitStateFunction);

% Compute the reliability index:
ReliabilityIndex = mean(LimitStateFunction)/std(LimitStateFunction);

%% Plot the results:
BinWidth = 2;

figure 
hold on
histogram(UltimateSystemStrength, "Normalization", "Probability", "BinWidth", BinWidth)
histogram(SystemLoad, "Normalization", "Probability", "BinWidth", BinWidth)
histogram(LimitStateFunction, "Normalization", "Probability", "BinWidth", BinWidth)
hold on
xlabel("$C_{System}$ or $Q_{System}$", "Interpreter", "Latex")
ylabel("PDF", "Interpreter", "Latex")
grid on
grid minor

% figure
% plot([0, SystemDisplacement', MaxNIncrements*DisplacementIncrement], [0, SystemStrength', SystemStrength(end)], ".-", "LineWidth", 1, "DisplayName", "System")
% hold on
% for i = 1:NSprings
%     plot([0, TPDisplacements{i}, MaxNIncrements*DisplacementIncrement], [0, TPStrengths{i}, TPStrengths{i}(end)], "LineWidth", 1, "DisplayName", strcat("Spring #", num2str(i)))
% end
% hold off
% title("$F$-$\Delta$ curves", "Interpreter", "Latex")
% xlabel("Displacement ($\Delta$)", "Interpreter", "Latex")
% ylabel("Force ($F$)", "Interpreter", "Latex")
% legend("Location", "Southeast")
% xlim([0, MaxNIncrements*DisplacementIncrement])
% grid on
% grid minor
