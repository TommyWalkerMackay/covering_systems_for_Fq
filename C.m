function [Cq] = C (q, U, Pascal)
% Usage: [Cq] = C (q, U, Pascal)

% Computes C(q), for integer q >= 2, as defined in the paper.
% If no third argument is passed, computes the value of Pascal using
% Pascal_triangle(max(q)+1).
% Uses BigFloat.
% Now vectorised.

% Checks for correct input.
if nargin < 2
    error('Not enough input arguments.')
end
if nargin < 3
    Pascal = [];
end

qmax = max(q, [], 'all');

if isempty(Pascal)
    Pascal = Pascal_triangle(qmax + 1);
end

if size(q, 1) ~= 1 || size(q, 2) ~= size(U, 2) || size(U, 1) <= qmax
    error('Incorrect dimensions for q and/or U.')
end
U = U(1:qmax+1, :);
if ~isnumeric(q)
    error('q must be numeric (not BigFloat).')
end

% Computes powers of (3*q-1)/(q-1)^2 as a BigFloat

qfun = (3 * q - 1) ./ BigFloat((q - 1) .^ 2);
qpowvec = powers_of(qfun, qmax, true, 'd');

Ureverse = BigFloat(zeros(size(U)));
for j = 1 : length(q)
    Ureverse(1 : q(j) + 1, j) = U(q(j) + 1 : -1 : 1, j);
end

Uqq = Ureverse(1, :);

qchoosew = Pascal(q+1, 1:size(U, 1));

Cq = sum(qchoosew' .* (Ureverse .* qpowvec) ./ Uqq, 1);

end