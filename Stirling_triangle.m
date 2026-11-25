function [S] = Stirling_triangle (sz)
% Usage: [S] = Stirling_triangle (sz)

% Computes a matrix of size sz x sz whose (n+1, k+1)-th element is S(n, k),
% where S(n, k) represents the Stirling number of the second kind.
% Uses the recurrence S(n+1, k) = k * S(n, k) + S(n, k-1).
% Uses BigFloat.

% Checks for correct input.
if ~isscalar(sz) || ~isreal(sz) || sz < 0 || floor(sz) ~= sz
    error('Input must be a non-negative integer.');
end

% Initialises S.

S = zeros(sz);
S(1, 1) = 1;
S = BigFloat(S);

% Applies the recurrence relation.
for rown = 2 : sz
    S(rown, :) = (0 : sz - 1) .* S(rown - 1, :) + ...
        [0, S(rown - 1, 1 : sz - 1)];
end

end