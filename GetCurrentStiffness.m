function [CurrentSpringStiffness] = GetCurrentStiffness(CurrentSpringElongation, SpringTPDisplacements, SpringTPStiffnessCoefficients, SpringInitialStiffness)
% Add the origin point to the turn-point displacements:
SpringTPDisplacements = [0, SpringTPDisplacements];

% Identify the linear interval number corresponding to the current elongation:
LocationIndex = find(CurrentSpringElongation >= SpringTPDisplacements, 1, "Last");

% Compute the stiffness in the current linear interval:
CurrentSpringStiffness = SpringTPStiffnessCoefficients(LocationIndex + 1)*SpringInitialStiffness;
end