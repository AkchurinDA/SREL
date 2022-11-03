% Author: Damir Akchurin
% Date: 11/1/22
% Version history:
% V1.0 - Damir Akchurin - First working version of the software (11/2/22)
% V1.1 - Damir Akchurin - Fixed an issue with finding parallel elements and added plotting function (11/3/22)

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
    % To simplify further comments, "element-containing node pairs" are refered to as "node pair".

    if size(ElemsInfo) ~= 1
        % Identify node pairs:
        [OriginalNodePairs, ~, GroupIndex] = unique(ElemsInfo(:, [2, 3]), 'Rows');
        NodePairs = OriginalNodePairs;
        NumElemsInEachNodePair = accumarray(GroupIndex, 1);

        % Identify node pairs that contain parallel elements which can be replaced with an equivalent
        % element:
        NewNodePairsStore = true;

        while ~isempty(NewNodePairsStore)
            % Preallocate:
            NumNodePairs = size(NodePairs, 1);
            AlreadyExistingNodePairsStore = cell(NumNodePairs, 1);
            NewNodePairsStore = cell(NumNodePairs, 1);

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
                    FormedNodePairs = zeros(NumCompatibleNodePairs, 2);

                    for j = 1:NumCompatibleNodePairs
                        FormedNodePairs(j, [1, 2]) = [CurrentNodePair(1), CompatibleNodePairs(j, 2)];
                    end
                else
                    FormedNodePairs = double.empty(0, 2);
                end

                % Categorize formed node pairs into new ones and the ones that already exist:
                AlreadyExistingNodePairs = intersect(NodePairs, FormedNodePairs, 'Rows');
                NewNodePairs = FormedNodePairs(~ismember(FormedNodePairs, AlreadyExistingNodePairs, 'Rows'), :);

                % Store new node pairs and the ones that already exist for the record:
                AlreadyExistingNodePairsStore{i} = AlreadyExistingNodePairs;
                NewNodePairsStore{i} = NewNodePairs;
            end

            % Add new node pairs to the existing element-containing node pairs:
            NewNodePairsStore = cell2mat(NewNodePairsStore);
            NumNewNodePairs = size(NewNodePairsStore, 1);
            NodePairs(1:end + NumNewNodePairs, :) = [NodePairs; NewNodePairsStore];
        end

        % Identify node pairs that contain parallel elements which can be replaced with an equivalent
        % element:
        AlreadyExistingNodePairsStore = cell2mat(AlreadyExistingNodePairsStore);
        if ~isempty(AlreadyExistingNodePairsStore)
            ParallelIndex = and(~ismember(OriginalNodePairs, AlreadyExistingNodePairsStore, 'Rows'), NumElemsInEachNodePair > 1);
            NodePairsWithParallelElems = OriginalNodePairs(ParallelIndex, :);
        else
            NodePairsWithParallelElems = OriginalNodePairs;
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
end

%% Plot the system:
figure
tiledlayout flow
for i = 1:numel(ElemsInfoStore)
    nexttile
    PlotHybridSpringSystem(ElemsInfoStore{i}, SysNodeI, SysNodeJ);
    title("Iteration " + (i - 1))
end
