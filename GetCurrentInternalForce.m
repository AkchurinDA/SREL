function [CurrentSpringAxialForce] = GetCurrentInternalForce(CurrentSpringElongation, SpringTPDisplacements, SpringTPStrengths, SpringTPStiffnessCoefficients, SpringInitialStiffness)
SpringTPDisplacements = [0, SpringTPDisplacements];
SpringTPStrengths = [0, SpringTPStrengths];

% Identify the linear interval number corresponding to the current elongation:
LocationIndex = find(CurrentSpringElongation >= SpringTPDisplacements, 1, "Last");

CurrentSpringAxialForce = SpringTPStrengths(LocationIndex) + SpringTPStiffnessCoefficients(LocationIndex + 1)*SpringInitialStiffness*(CurrentSpringElongation - SpringTPDisplacements(LocationIndex));
end