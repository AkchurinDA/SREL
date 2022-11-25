function PlotForceDisplacementCurves(SpringInflectionDisplacements, SpringInflectionStrength, SystemDisplacementRecord, SystemSupportedLoadRecord, NumSprings, NumIterations, DisplacementIncrement)
figure("Units", "Inches", "Position", [1, 1, 18, 6])

subplot(1, 3, 1)
hold on
for i = 1:NumSprings
    plot([0, SpringInflectionDisplacements{i}, NumIterations*DisplacementIncrement], [0, SpringInflectionStrength{i}, SpringInflectionStrength{i}(end)], "LineWidth", 1, "DisplayName", strcat("Spring #", num2str(i)))
end
hold off
title("$F$-$\Delta$ curves for individual springs", "Interpreter", "Latex")
xlabel("Displacement ($\Delta_{Spring}$)", "Interpreter", "Latex")
ylabel("Force ($F_{Spring}$)", "Interpreter", "Latex")
legend("Location", "Best")
grid on
grid minor

subplot(1, 3, 2)
plot(SystemDisplacementRecord, SystemSupportedLoadRecord, 'k-', "LineWidth", 1)
title("$F$-$\Delta$ curve for the system", "Interpreter", "Latex")
xlabel("Displacement ($\Delta_{System}$)", "Interpreter", "Latex")
ylabel("Force ($F_{System}$)", "Interpreter", "Latex")
grid on
grid minor

subplot(1, 3, 3)
hold on
for i = 1:NumSprings
    plot([0, SpringInflectionDisplacements{i}, NumIterations*DisplacementIncrement], [0, SpringInflectionStrength{i}, SpringInflectionStrength{i}(end)], "LineWidth", 1, "DisplayName", strcat("Spring #", num2str(i)))
end
plot(SystemDisplacementRecord, SystemSupportedLoadRecord, 'k-', "LineWidth", 1, "DisplayName", "System")
hold off
title("Comparison of $F$-$\Delta$ curves", "Interpreter", "Latex")
xlabel("Displacement ($\Delta$)", "Interpreter", "Latex")
ylabel("Force ($F$)", "Interpreter", "Latex")
legend("Location", "Best")
grid on
grid minor