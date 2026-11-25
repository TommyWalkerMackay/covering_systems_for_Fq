function [B_tilde] = B_tilde_vector (n_max)
% Usage: [B_tilde] = B_tilde_vector (n_max)

% Computes a column vector of length (n_max + 1) whose (n+1)-th element is
% B_tilde(n), where B_tilde(n) is the n-th complementary Bell number.
% Uses BigFloat.

% Checks for correct input.
if ~isscalar(n_max) || ~isreal(n_max) || n_max < 0 || floor(n_max) ~= n_max
    error('Input must be a non-negative integer.');
end

% Computes B_tilde using the definition of B_tilde in terms of Stirling 
% numbers of the second kind.
S = Stirling_triangle(n_max + 1);
B_tilde = S * (-1) .^ (0 : n_max)';

end