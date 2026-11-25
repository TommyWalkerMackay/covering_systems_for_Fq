function [V] = V_bound (q, delta_vec, C)
% Usage: [V] = V_bound (q, delta_vec, C)

% Computes an upper bound on V using the method described in the proof of
% Theorem 1.3. delta_vec should begin with delta(2) and continue until
% delta(10), or whatever maximum value of d we wish to use in the explicit
% computation as opposed to the cruder bound.

% Checks for correct input.

if ~isvector(delta_vec)
    error('delta_vec must be a vector.')
end
if size(delta_vec, 1) > 1
    delta_vec = delta_vec';
end
if any(floor(q) ~= q, 'all')
    error('q must be an integer.')
end

szq = size(q);
if ndims(C) ~= ndims(q) || any(size(C) ~= szq)
    error('Sizes of C and q must match.')
end
q = q(:);
C = C(:);

d_max = length(delta_vec) + 1;

% Computes N_q(d) for 2 <= d <= d_max
Nq = Nq_vector(q, d_max);

% Computes g(d) for 2 <= d <= d_max

qpowvec = powers_of(q, max(d_max, 4), false, 'r');
qm1powvec = powers_of(q - 1, 6, false, 'r');

g = (3 * qpowvec(:, 2:d_max) - 1) ./ (qpowvec(:, 2:d_max) - 1) ./ (qpowvec(:, 2:d_max) - 1);

% Computes (1 + (1 - delta(d))^-1 g(d)) ^ N_q(d) for 2 <= d <= d_max.
% Uses oneplusxton to get a rational upper bound on this quantity.

keyfactor = oneplusxton(g ./ (1 - delta_vec), Nq);

VdoverC = 1 ./ (4 * delta_vec .* (3 * qpowvec(:, 2 : d_max) - 1)) .* ...
    (keyfactor - 1) .* cumprod([ones(length(q), 1), keyfactor(:, 1 : end - 1)], 2);

remainingVterm = q .* (qpowvec(:, 4) + 26*qpowvec(:, 3) + 66*qpowvec(:, 2) + 26*q + 1) ./ qm1powvec(:, 6) - ...
    sum((1 : d_max) .^ 5 ./ qpowvec(:, 1:d_max), 2);

V = C .* (sum(VdoverC, 2) + remainingVterm);

V = reshape(V, szq);

end