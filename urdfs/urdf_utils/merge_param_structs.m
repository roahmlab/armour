function [out] = merge_param_structs(x,y)
    % given two structs x and y, merge and return the ouput.
    out = cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
end

