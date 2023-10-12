function [newzono] = zonotope_slice(zono, slice_dim, slice_pt)
% function for taking a slice of a zonotope along the corresponding
% dimensions...
% for now, only works if there's one generator in each slice_dim

if isempty(slice_dim) && isempty(slice_pt)
    newzono = zono;
    return;
elseif isempty(slice_dim) || isempty(slice_pt)
    warning('either slice_dim or slice_value is nonempty... leave both empty if not slicing. returning original zono');
    newzono = zono;
    return;
end

if (size(slice_dim, 1) ~= size(slice_pt, 1))
  error('number of rows of slice_dim and slice_pt are unequal');
end

if ~ (size(slice_pt, 2) == 1 || size(slice_pt, 2) == 2)
    error('Slice point should be a column vector or two column vectors specifying an interval');
end

Z = zono.Z;
c = Z(:, 1);
G = Z(:, 2:end);

slice_idx = [];
for i = 1:length(slice_dim)
   myidxs = find(G(slice_dim(i), :) ~= 0);
   if length(myidxs) ~= 1
       if length(myidxs) == 0
           error('No generator for slice index');
       else
           error('More than one generator for slice index');
       end
   end
   slice_idx(i, 1) = myidxs;
end

slice_c = c(slice_dim, 1);
slice_G = G(slice_dim, slice_idx);
for i = 1:size(slice_pt, 2)
  slice_lambda(:, i) = slice_G\(slice_pt(:, i) - slice_c);
  if any(abs(slice_lambda(:, i)) > 1)
    error('Slice point is outside bounds of reach set, and therefore is not verified');
  end
end

newc = c;
newG = G;
if size(slice_lambda, 2) == 1
  newG(:, slice_idx) = [];
  newc = newc + G(:, slice_idx)*slice_lambda;
elseif size(slice_lambda, 2) == 2
    for i = 1:size(slice_dim, 1)
      tmp_lb = newG(:, slice_idx(i))*slice_lambda(i, 1);
      tmp_ub = newG(:, slice_idx(i))*slice_lambda(i, 2);

      newc = newc + (tmp_lb + tmp_ub)/2; % shift center 
      newG(:, slice_idx(i)) = (tmp_ub - tmp_lb)/2; % scale generator
    end
    newG(:, ~any(newG)) = []; % remove zero columns (i.e. those caused...
    % by slicing by an interval whose lower bound equals its upper bound)
end

newzono = zonotope([newc, newG]);

end

% % function for taking a slice of a zonotope along the corresponding
% % dimensions...
% % for now, only works if there's one generator in each slice_dim
% 
% % if size(slice_pt, 2) ~= 1
% %     error('Slice point should be a column vector');
% % end
% 
% if (size(slice_dim, 1) ~= size(slice_pt, 1))
%   error('number of rows of slice_dim and slice_pt are unequal');
% end
% 
% if ~ (size(slice_pt, 2) == 1 || size(slice_pt, 2) == 2)
%     error('Slice point should be a column vector or two column vectors specifying an interval');
% end
% 
% Z = get(zono, 'Z');
% c = Z(:, 1);
% G = Z(:, 2:end);
% 
% slice_idx = [];
% for i = 1:length(slice_dim)
%    myidxs = find(G(slice_dim(i), :) ~= 0);
%    if length(myidxs) ~= 1
%        if length(myidxs) == 0
%            error('No generator for slice index');
%        else
%            error('More than one generator for slice index');
%        end
%    end
%    slice_idx(i, 1) = myidxs;
% end
% 
% slice_c = c(slice_dim, 1);
% slice_G = G(slice_dim, slice_idx);
% slice_lambda = slice_G\(slice_pt - slice_c);
% if size(slice_lambda, 2) > 1
%     error('slice_lambda is not 1D');
% end
% if any(abs(slice_lambda) > 1)
%     error('Slice point is outside bounds of reach set, and therefore is not verified');
% end
% 
% newG = G;
% newG(:, slice_idx) = [];
% newc = c + G(:, slice_idx)*slice_lambda;
% 
% newzono = zonotope([newc, newG]);
% 
% end

