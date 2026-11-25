function [U] = U_vector (q, B_tilde, Pascal)
% Usage: [U] = U_vector (q, B_tilde, Pascal)

% Computes a column vector of length q+1 whose (s+1)-th element is U(q, s),
% where U(q, s) is the uncovered space function (taking nu -> infinity).
% If no second argument is passed, computes the value of B_tilde using
% B_tilde_vector(q).
% If no third argument is passed, computes the value of Pascal using
% Pascal_triangle(q+1).
% Uses BigFloat.
% Now adapted for vector q input.

% Checks for correct input.
if nargin < 2
    B_tilde = [];
end
if nargin < 3
    Pascal = [];
end
if ~isvector(q) || ~isnumeric(q)
    error('q must be a numeric vector (not a BigFloat).')
end
qmax = max(q, [], 'all');
if isempty(Pascal)
    Pascal = Pascal_triangle(qmax + 1);
end
if isempty(B_tilde)
    B_tilde = B_tilde_vector(qmax);
end
if ~isvector(B_tilde) || length(B_tilde) < qmax + 1
    error('Second argument must be a vector of length at least q+1.')
end

q = q(:)';

% Truncates B_tilde, Pascal to the required size.
B_tilde = B_tilde(1 : qmax+1);
Pascal = Pascal(1 : qmax+1, 1 : qmax+1);

% Computes powers of 1/(q-1) as a BigFloat
tic
disp('Computing powers of 1/(q-1)...')
qpowvec = powers_of(1 ./ BigFloat(q - 1), qmax, true, 'd');
toc

% Computes U.
tic
disp('Computing matrix product...')
x = B_tilde .* qpowvec;
U = Pascal * x;
toc

end