function [tf] = is_discrete_cov_sys (coords, dims)
% Usage: [tf] = is_discrete_cov_sys (coords, dims)

% Takes cell input coords and vector input dims.
% Each row of coords represents a hyperplane, with a fixed co-ordinate
% represented by that co-ordinate and a varying coordinate by ':'.
% Returns true if the hyperplanes represented by the rows of coords are a
% distinct covering system for [dims(1)] x [dims(2)] x ... x [dims(end)].

if isscalar(dims)
    dims = [dims, 1];
end

covered = false(dims);

for n = 1 : size(coords, 1)
    hyp = coords(n, :);
    covered(hyp{:}) = true;
end

covers = all(covered, 'all');

fixed_coords = true(size(coords));
for j = 1 : numel(coords)
    if ischar(coords{j})
        fixed_coords(j) = false;
    end
end

% Assigns a unique number to each possible fixed co-ordinate set.
fixreps = fixed_coords * 2 .^ (1 : size(fixed_coords, 2))';


% False if any two hyperplanes have the same fixed co-ordinates.
norepeatdims = ~any(fixreps == fixreps' & ~eye(length(fixreps)), 'all');
% False if any hyperplane is trivial.
notrivialhyps = ~any(fixreps == 0, 'all');

if covers && notrivialhyps && norepeatdims
    tf = true;
else
    tf = false;
end
end
