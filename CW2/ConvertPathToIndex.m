% Generate image sequence from given paths (no duplication between 2 frames)
function indexSequence = ConvertPathToIndex(paths)
disp("@Convert Path to Index");
tic  
    indexSequence = [];
    lastIndex = 0;
    for p = 1:length(paths)
        path = paths{p};
        for i = 1:length(path)
            if (path(i) ~= lastIndex)
                indexSequence = [indexSequence path(i)];
                lastIndex = path(i);
            end
        end
    end
toc
end