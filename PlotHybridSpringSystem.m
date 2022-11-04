% Author: Damir Akchurin
% Date: 11/3/22

function PlotHybridSpringSystem(ElemsInfo, SysNodeI, SysNodeJ)
% Create a graph using nodes of each element:
G = graph(ElemsInfo(:, 2), ElemsInfo(:, 3));

% Remove 0-degree nodes from the graph:
G = rmnode(G, find(degree(G) == 0));

% Set parameters of the graph:
ElementWidth = 1.5;
ElementColor = 'K';
NodeSize = 10;
NodeColor = '#0072BD';
SysNodeMarker = '^';
SysNodeColor = 'R';

% Define node names:
NodeLabelsNumerical = unique(ElemsInfo(:, [2, 3]));
NodeLabelsCell = cell(numel(NodeLabelsNumerical), 1);
for i = 1:numel(NodeLabelsNumerical)
    NodeLabelsCell{i} = num2str(NodeLabelsNumerical(i));
end

% Assign node names:
G.Nodes.Name = NodeLabelsCell;

% Define edge labels:
ElemsInfo = sortrows(ElemsInfo, [2, 3]);
EgdeLabels = strcat(num2str(ElemsInfo(:, 1)), repmat(", $\beta$ = ", size(ElemsInfo, 1), 1), num2str(ElemsInfo(:, 4)));

% Plot the graph:
H = plot(G, 'NodeColor', NodeColor, 'MarkerSize', NodeSize, 'EdgeColor', ElementColor, 'LineWidth', ElementWidth,...
    'EdgeLabel', EgdeLabels, 'Interpreter', 'Latex',...
    'Layout', 'SubSpace');

% Highlight the system nodes:
SysNodeIID = findnode(G, num2str(SysNodeI));
SysNodeJID = findnode(G, num2str(SysNodeJ));
highlight(H, [SysNodeIID, SysNodeJID], 'Marker', SysNodeMarker, 'NodeColor', SysNodeColor)

% Add additional info:
title('Graphical representation of the system')
end