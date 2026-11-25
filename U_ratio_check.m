function [success, minratios, ratios_error, argminratios] = U_ratio_check (q, U)
% Usage: [success, minratios, ratios_error, argminratios] = U_ratio_check (q, U)

% Inputs q and column vector U whose (s+1)-th element is U(q, s).
% Returns true if for all 1 <= s <= q, U(q, s) / U(q, s-1) > 1/q, and
% returns false otherwise.
% Now also outputs:
% minratios = min_{s leq q} U(q,s) / U(q,s-1),
% ratios_error = the maximum difference between upper and lower bounds on
% the ratios
% argminratios = argmin_{s leq q} U(q,s) / U(q,s-1)

if ~isvector(q) || ~isnumeric(q)
    error('q must be a numeric vector.')
end
q = q(:)';

sz = size(U, 1);
if sz <= max(q, [], 'all')
    error('Size of U must be greater than the maximum value in q.');
end

upperU = U.upper();
lowerU = U.lower();
zeroinU = containsZero(U);

ratios = lowerU(2:sz, :) ./ upperU(1:sz-1, :);
upper_ratios = upperU(2:sz, :) ./ lowerU(1:sz-1,:);
outofbounds = q < (1:sz-1)';
low_ratio = ~ (ratios > 1 ./ q) | zeroinU(1:sz-1, :) | isinf(upperU(1:sz-1, :));
success = ~any(low_ratio & ~outofbounds, 'all');

ratios(outofbounds) = Inf;
minratios = min(ratios, [], 1);
ratios_error = max(upper_ratios(~outofbounds) - ratios(~outofbounds), [], 'all');

argminratios_log = ratios == minratios;
argminratios = zeros(size(q));
for j = 1 : length(q)
    argminratios(j) = find(argminratios_log(:, j), 1);
end

end