function [a] = oneplusxton (x, n)
% Usage: [a] = oneplusxton (x, n)

% Outputs (1 + x) .^ n, where x and n are BigFloats.

if any(double(x) <= -1, 'all')
    error('x must be strictly greater than -1.')
end

x = BigFloat(x);
n = BigFloat(n);

max_dim = max(ndims(x), ndims(n));
dims = max(size(x, 1:max_dim), size(n, 1:max_dim));
nelts = prod(dims);
x = x .* ones(dims);
n = n .* ones(dims);

lo_x = BigFloat(x.lo_m(), x.lo_e());
hi_x = BigFloat(x.hi_m(), x.hi_e());
lo_n = BigFloat(n.lo_m(), n.lo_e());
hi_n = BigFloat(n.hi_m(), n.hi_e());

sqrt2 = BigFloat(sqrt(2)) + BigFloat(0);

a = BigFloat(zeros(size(x)));

for k = 1:4
    if k == 1 || k == 2
        this_x = lo_x;
    else
        this_x = hi_x;
    end
    if k == 1 || k == 3
        this_n = lo_n;
    else
        this_n = hi_n;
    end

    oneplusx = 1 + this_x;
    oneplusx_m = oneplusx.lo_m();
    oneplusx_e = double(oneplusx.lo_e());

    % We shall factor (1+x) = (1+r) .* sqrt(2) .^ (n_sqrt2_factors).
    n_sqrt2_factors = 2 * oneplusx_e - (oneplusx_m < sqrt(2)/2) - (oneplusx_e >= 1);
    sqrt2s_and_ones = BigFloat(ones(size(oneplusx)));
    sqrt2s_and_ones(mod(n_sqrt2_factors, 2) == 1) = sqrt2;
    twopow = ceil(n_sqrt2_factors ./ 2);

    % Define r, nm.
    r = oneplusx .* sqrt2s_and_ones ./ ...
        BigFloat(ones(size(oneplusx)), twopow) - 1;
    r(n_sqrt2_factors == 0) = this_x(n_sqrt2_factors == 0);

    % Ensures abs(min(x, 1/2) .* n_red) <= 1/2
    n_exp = max(this_n.lo_e() + min(this_x.lo_e(), 0), 0) + 1;
    n_red = BigFloat(this_n.lo_m(), this_n.lo_e() - n_exp);

    % Initialize matrix to take values of (1+x) .^ nm
    b = BigFloat(ones(size(r)));

    for j = 1:2
        if j == 1
            this_r = r;
        else
            this_r = sqrt2 - 1;
        end

        % Compute binomial expansion of (1 + this_r) .^ nm.
        nterms = 60;
        summands = BigFloat(zeros([dims, nterms]));
        summands(1:nelts) = BigFloat(ones(nelts, 1));
        % Calculate the binomial expansion terms
        for term = 2:nterms
            idx = cell(1, max_dim+1);
            idx(:) = {':'};
            idx{end} = term;
            idx2 = cell(1, max_dim+1);
            idx2(:) = {':'};
            idx2{end} = term-1;
            summands(idx{:}) = summands(idx2{:}) .* (n_red + 2 - term) ...
                ./ (term - 1) .* this_r;
        end
        % This is an upper bound on the leftover sum provided that
        % abs(r .* n_red) <= 1/2. This is guaranteed by the definition of
        % n_red in terms of x and the fact that abs(r) <= abs(x) (unless
        % n_sqrt_2_factors < 0 in which case this won't be counted).
        max_remainder = BigFloat(1, 1 - nterms);

        binom = sum(summands, max_dim + 1) + ...
            combine_range(max_remainder, - max_remainder);
        if j == 1
            b = b .* binom;
        else
            for it = 1 : max(abs(n_sqrt2_factors), [], 'all')
                b(n_sqrt2_factors >= it) = b(n_sqrt2_factors >= it) .* ...
                    binom(n_sqrt2_factors >= it);
                b(n_sqrt2_factors <= -it) = b(n_sqrt2_factors <= -it) ./ ...
                    binom(n_sqrt2_factors <= -it);
            end
        end
    end

    this_a = b;
    for it = 1 : max(n_exp, [], 'all')
        to_square = n_exp >= it;
        this_a(to_square) = this_a(to_square) .* this_a(to_square);
    end

    if k == 1
        a = this_a;
    else
        a = combine_range(a, this_a);
    end
    
end
end
