function [pows] = powers_of (v, max_pow, include_0, rd)
% Usage: [pows] = powers_of (v, max_pow, include_0, rd)

% Computes pointwise powers of vector v.
% max_pow = maximum power to compute
% Includes v .^ 0 = 1 if include_0 is set to true, doesn't otherwise.
% If rd = 'd' then the powers are arranged vertically, if rd = 'r' then
% they are arranged horizontally.

if nargin < 3
    include_0 = [];
end
if nargin < 4
    rd = [];
end

if isempty(include_0)
    include_0 = false;
end

% Determines whether we are going right or down.
if isempty(rd)
    if size(v, 2) > 1
        rd = 'd';
    else
        rd = 'r';
    end
end
if ~isvector(v)
    error('v must be a vector.')
end

v = BigFloat(v);

npows = max_pow + include_0;
sz = size(v);
if rd == 'd'
    sz(1) = npows;
else
    sz(2) = npows;
end

pows = BigFloat(zeros(sz));

if ~include_0
    if rd == 'd'
        pows(1, :) = v; % Assign v for the 1st power
    else
        pows(:, 1) = v; % Assign v for the 1st power
    end
else
    if rd == 'd'
        pows(1, :) = BigFloat(1); % Assign 1 for the 0th power
    else
        pows(:, 1) = BigFloat(1); % Assign 1 for the 0th power
    end
end

for i = 2:npows
    if rd == 'd'
        pows(i, :) = v .* pows(i-1, :);
    else
        pows(:, i) = v .* pows(:, i-1);
    end
end
end