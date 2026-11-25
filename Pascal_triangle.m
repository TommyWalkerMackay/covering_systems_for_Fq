function [P] = Pascal_triangle (sz)
% Usage: [P] = Pascal_triangle (sz)

% Computes a matrix of size sz x sz whose (n+1, k+1)-th element is 
% n choose k.
% Uses the recurrence n+1 choose k = n choose k + n choose k-1.
% Uses BigFloat.

% Checks for correct input.
if ~isscalar(sz) || ~isreal(sz) || sz < 0 || floor(sz) ~= sz
    error('Input must be a non-negative integer.');
end
sz = double(sz);

% Initialize P.
P = zeros(sz);
P(1, 1) = 1;
P = BigFloat(P);

% Applies the recurrence relation.
for rown = 2 : sz
    P(rown, :) = P(rown - 1, :) + [0, P(rown - 1, 1 : sz - 1)];
end

% Sets all elements that should be 0 exactly to 0.
P((1:sz)' < 1:sz) = BigFloat(0);

end