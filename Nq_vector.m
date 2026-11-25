function [Nq] = Nq_vector (q, d_max)
% Usage: [Nq] = Nq_vector (q, d_max)

% Outputs a row vector of length d_max - 1 whose (d - 1)-th element is
% N_q(d).

% Checks for correct input.
if ~isvector(q)
    error('q must be a vector.')
end
if ~isscalar(d_max)
    error('d_max must be a scalar.')
end
if d_max < 2
    error('Maximum d must be at least 2.')
end
if any(floor(q) ~= q, 'all')
    error('q must be an integer.')
end

q = q(:);

% Initialize Nq.
Nq = zeros(length(q), d_max);
Nq(:, 1) = q;
Nq = BigFloat(Nq);

% Computes powers of q.

qpowvec = powers_of(q, d_max, false, 'r');

% Calculate N_q(d) for each d from 2 to d_max.
for d = 2 : d_max
    Nqd = (qpowvec(:, d) - Nq .* (1 : d_max) * double(mod(d, (1 : d_max)') == 0)) ./ d;
    Nq(:, d) = Nqd;
end

% Delete first element of Nq.
Nq = Nq(:, 2 : end);

end