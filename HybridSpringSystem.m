% Author: Damir Akchurin
% Date: 1/1/22
% Version history:
% V1.0 - First working version of the software (1/2/22)

clc
close all
clear variables

%% Define element information array:
% [Element # | Node I # | Node J # | Reliability index]
ElemsInfo = [
    1, 1, 4, 2.5;
    2, 4, 3, 2.5;
    3, 1, 2, 2.5;
    4, 1, 2, 2.5;
    5, 2, 3, 2.5;
    6, 1, 3, 2.5;
    7, 3, 5, 2.5];

%% Define system nodes:
SysNodeI = 1;
SysNodeJ = 5;

%% Calculate system reliability index:
% Store element information array to track the results:
Counter = 1;
ElemsInfoStore{Counter} = ElemsInfo;

% Sort node numbering to simplify further calculations:
ElemsInfo(:, [2, 3]) = sort(ElemsInfo(:, [2, 3]), 2);

while size(ElemsInfo, 1) ~= 1
    %% Check if any elements are in series:
    % NOTE: 
    % - If a node has two elements connected to it and if it is not a system node, then these
    % elements are in series and can be replaced with an equivalent element
    % - If a node has 2 elements connected to it and if it is a system node, then these elements are
    % not in series, but in parallel.

    % Calculate the number of elements connected to each node:
    [NumConnectedElems, Nodes] = groupcounts(reshape(ElemsInfo(:, [2, 3]), [], 1));

    % Identify to what nodes have elements connected in series:
    SeriesIndex = NumConnectedElems == 2;
    SeriesNodes = Nodes(SeriesIndex);

    if ~isempty(SeriesNodes)
        for i = 1:numel(SeriesNodes)
            % Make sure that the node is not a system node, beucase it might indicate that the node
            % might connect two elements in parallel
            if i ~= SysNodeI || i ~= SysNodeJ 
                % Identify two elements that are in series:
                [SeriesElemsIndexR, SeriesElemsIndexC] = find(ElemsInfo(:, [2, 3]) == SeriesNodes(i));

                % Calculate equivalent probability of failure:
                FailureProbability = cdf("Normal", -ElemsInfo(SeriesElemsIndexR, 4), 0, 1);
                EqFailureProbability = 1 - prod(1 - FailureProbability);

                % Calculate equivalent reliability index:
                EqReliabilityIndex = -icdf("Normal", EqFailureProbability, 0, 1);

                % Update the element information array:
                NewElemNum = ElemsInfo(SeriesElemsIndexR(1), 1);
                NewNodeI = ElemsInfo(SeriesElemsIndexR(1), 1 + (3 - SeriesElemsIndexC(1)));
                NewNodeJ = ElemsInfo(SeriesElemsIndexR(2), 1 + (3 - SeriesElemsIndexC(2)));
                NewNodes = sort([NewNodeI, NewNodeJ]);
                ElemsInfo(SeriesElemsIndexR, :) = [];
                ElemsInfo(size(ElemsInfo, 1) + 1, :) = [NewElemNum, NewNodes, EqReliabilityIndex];
            end
        end
    end
    
    Counter = Counter + 1;
    ElemsInfoStore{Counter} = ElemsInfo;

    %% Check if any elements are in parallel:
    % NOTE: If combinations of element-containing node pairs can form other element-containing node
    % pairs, then the former node pairs, given that they contain more than one element, must contain
    % elements that are connected in parallel and can be replaced with an equivalent element.

    % Identify element-containing node pairs:
    [NodePairs, ~, GroupIndex] = unique(ElemsInfo(:, [2, 3]), 'Rows');
    NumElemsInEachNodePair = accumarray(GroupIndex, 1);

    % Identify node pairs that contain parallel elements which can be replaced with an equivalent
    % element:
    NumNodePairs = size(NodePairs, 1);
    BannedNodePairs = cell(NumNodePairs, 1);
    for i = 1:NumNodePairs
        % Extract a node pair:
        CurrentNodePair = NodePairs(i, :);

        % Find compatible node pairs:
        CompatibleNodePairsIndex = CurrentNodePair(2) == NodePairs(:, 1);
        CompatibleNodePairs = NodePairs(CompatibleNodePairsIndex, :);

        % Find nodes pairs that can be formed from combining the current node pair and node pairs
        % that are compatible with it:
        if ~isempty(CompatibleNodePairs)
            NumCompatibleNodePairs = size(CompatibleNodePairs, 1);
            NewlyFormedNodePairs = zeros(NumCompatibleNodePairs, 2);

            for j = 1:NumCompatibleNodePairs
                NewlyFormedNodePairs(j, [1, 2]) = [CurrentNodePair(1), CompatibleNodePairs(j, 2)];
            end
        else
            NewlyFormedNodePairs = [];
        end
        
        % Add newly formed node pairs to a ban list
        BannedNodePairs{i} = NewlyFormedNodePairs;
    end
    
    BannedNodePairs = cell2mat(BannedNodePairs);
    if ~isempty(BannedNodePairs)
        ParallelIndex = and(~ismember(NodePairs, BannedNodePairs, 'Rows'), NumElemsInEachNodePair > 1);
        NodePairsWithParallelElems = NodePairs(ParallelIndex, :);
    else 
        NodePairsWithParallelElems = NodePairs;
    end

    if ~isempty(NodePairsWithParallelElems)
        for i = 1:size(NodePairsWithParallelElems, 1)
            % Identify elements that are in parallel:
            ParallelElemsIndex = ismember(ElemsInfo(:, [2, 3]), NodePairsWithParallelElems(i, :), 'Rows');

            % Calculate equivalent probability of failure:
            FailureProbability = cdf("Normal", -ElemsInfo(ParallelElemsIndex, 4), 0, 1);
            EqFailureProbability = prod(FailureProbability);

            % Calculate equivalent reliability index:
            EqReliabilityIndex = -icdf("Normal", EqFailureProbability, 0, 1);

            % Update the element information array:
            ParallelElems = ElemsInfo(ParallelElemsIndex, 1);
            NewElemNum = ParallelElems(1);
            ElemsInfo(ParallelElemsIndex, :) = [];
            ElemsInfo(size(ElemsInfo, 1) + 1, :) = [NewElemNum, NodePairsWithParallelElems(i, :), EqReliabilityIndex];
        end
    end

    Counter = Counter + 1;
    ElemsInfoStore{Counter} = ElemsInfo;
end

