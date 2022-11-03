% Author: Damir Akchurin
% Date: 11/3/22

function PlotHybridSpringSystem(ElemsInfo, SysNodeI, SysNodeJ)
% Create a graph using nodes of each element:
G = graph(ElemsInfo(:, 2), ElemsInfo(:, 3));

% Set parameters of the graph:
ElementWidth = 1.5;
ElementColor = 'K';
NodeSize = 10;
NodeColor = '#0072BD';
SysNodeMarker = '^';
SysNodeColor = 'R';

% Plot the graph:
H = plot(G, 'NodeColor', NodeColor, 'MarkerSize', NodeSize, 'EdgeColor', ElementColor, 'LineWidth', ElementWidth,...
    'Layout', 'SubSpace', 'EdgeLabel', strcat(num2str(ElemsInfo(:, 1)), repmat(", $\beta$ = ", size(ElemsInfo, 1), 1), num2str(ElemsInfo(:, 4))), 'Interpreter', 'Latex');

% Highlight the system nodes:
highlight(H, [SysNodeI, SysNodeJ], 'Marker', SysNodeMarker, 'NodeColor', SysNodeColor)

% Add additional info:
title('Graphical representation of the system')
end