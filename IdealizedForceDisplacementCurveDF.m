% Author: Damir Akchurin
% Date: 11/7/22

function CurrentForce = IdealizedForceDisplacementCurveDF(CurrentDisplacement, Stiffness, Coefficients, TurnPointForces)
% PURPOSE:
% To model the idealized force-displacement
% INPUT:
% CurrentDisplacement - diplacement of an element at a given instant
% Stiffness - initial stiffness (slope) of the first linear interval of the force-displacement curve
% Coefficients - coefficients that are used to obtain stiffnesses (slope) at each linear interval of the
% force-displacement curve
% TurnPointForces - forces at the intersections between each linear interval
% OUTPUT:
% CurrentForce - force within the element at a given instant

% Calculate displacements at turnover points:
TurnPointDisplacements = zeros(1, numel(Coefficients));
for i = 1:numel(Coefficients)
    if i == 1
        TurnPointDisplacements(i) = TurnPointForces(i)/(Coefficients(i)*Stiffness);
    else
        TurnPointDisplacements(i) = TurnPointDisplacements(i - 1) + (TurnPointForces(i) - TurnPointForces(i - 1))/(Coefficients(i)*Stiffness);
    end
end

% Add the origin as the very first turnover point:
TurnPointForces = [0, TurnPointForces];
TurnPointDisplacements = [0, TurnPointDisplacements];

% Find under which interval of the force-displacement curve the current displacement falls:
IntervalIndex = nnz(CurrentDisplacement > TurnPointDisplacements) + 1;

% Calculate the current force:
CurrentForce = TurnPointForces(IntervalIndex - 1) + Coefficients(IntervalIndex)*Stiffness*(CurrentDisplacement - TurnPointDisplacements(IntervalIndex - 1)); 
end