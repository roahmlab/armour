function [id, id_names, id_add] = update_id(id, id_names, n, label, varargin)
% want to keep track of all generator coefficients we have seen, and label
% each with a different name.
    id_len = length(id);
    id_add = (id_len+1:id_len+n)';
    id = [id; id_add];
    for i = 1:n
        id_names = [id_names; sprintf('%s_%03d', label, i)];
%         if isempty(varargin)
%             id_names = [id_names; sprintf('%s_%03d_%03d', label, i, 0)];
%         else
%             id_names = [id_names; sprintf('%s_%03d_%03d', label, i, varargin{1})];
%         end
    end
end

