function MAIN (latex)
% Usage: MAIN (latex)

% Performs all of the calculations required for the proof.

starttime = datetime;

if nargin < 1
    output_latex = false; % If set to true prints output in a way that can
                          % easily be copied to a LaTeX file.
else
    output_latex = latex;
end

% Checks the covers in Table 1 are correct.

disp('Checking the covering systems given in Table 1...')
check_cov_systems(output_latex);

% Computes miscellaneous approximations.

fprintf('\n0.505 e (e^1.01 sinh e^1.01 + e^2.02 cosh e^1.01) = ')
fprintf('%.8g.\n', 0.505 * exp(1) * (exp(1.01) * sinh(exp(1.01)) + exp(2.02) * cosh(exp(1.01))))
fprintf('e^(1 - e^1.01) = ')
fprintf('%.8g.\n', exp(1 - exp(1.01)))
fprintf('80 (1 + 3/700 + 2/700^2) e^3 = ')
fprintf('%.8g.\n', 80 * (1 + 3/700 + 2/700^2) * exp(3))
fprintf('1620 \\sum_{d=2}^{Inf} d^5/250^d = ')
f = @(x)x.*(x.^4+26*x.^3+66*x.^2+26*x+1)./(x-1).^6;
fprintf('%.8g.\n', 1620 * (f(250) - 1/250))

% The rest of this program deals with computing key values.

switch 1
    case 1 % Compute U and C.

    qmin = 7;
    qmax = 700;
    
    % Initialize q values
    
    qvalues = [];
    for q = qmin : qmax
        qfacts = factor(q);
        if all(qfacts == qfacts(1))
            % Then q is a prime power
            qvalues(end + 1) = q;
        end
    end
    
    % Compute B_tilde and Pascal
    tic
    fprintf('\nComputing B_tilde...\n')
    B_tilde = B_tilde_vector(qmax);
    toc
    tic
    fprintf('\nComputing Pascal_triangle...\n')
    Pascal = Pascal_triangle(qmax + 1);
    toc
    
    % Perform calculations for all q simultaneously.

    fprintf('\nComputing U...\n')
    U = U_vector(qvalues, B_tilde, Pascal);
    fprintf('\nChecking ratios U(q, s) / U(q, s-1)...\n')
    [success, minratios, ratios_error, argminratios] = U_ratio_check(qvalues, U);
    if success
        disp('Ratios are acceptable.');
    else
        disp('Ratios are incorrect.');
        keyboard;
    end
    tic
    fprintf('\nComputing C(q)...\n')
    C_vector = C(qvalues, U, Pascal);
    toc
    fprintf('Checking C(q)...\n')
    if all(C_vector.upper() < 1620, 'all')
        disp('C(q) < 1620 for all q.')
    else
        disp('C(q) may be >= 1620 for some q.')
        keyboard;
    end

    case 2 % Import U and C from a previous run.
    str = importdata("Covering_systems_data.mat");
    C_vector = str.C_vector;
    qvalues = str.qvalues;
    minratios = str.minratios;
    ratios_error = str.ratios_error;
    argminratios = str.argminratios;
end

newqmax = 250;
newqvalues = qvalues(qvalues <= newqmax);

tic
fprintf('\nComputing V bounds...\n')
delta_vec = [1/BigFloat(4), 1/BigFloat(2) .* ones(1, 8)];
thisC = C_vector(qvalues <= newqmax);
V_bounds = V_bound(newqvalues, delta_vec, thisC);
toc

V_upper = V_bounds.upper();
big_V_upper = NaN(size(qvalues));
big_V_upper(1 : length(V_upper)) = V_upper;

if all(argminratios == qvalues, 'all')
    fprintf('\nU(q,s)/U(q,s-1) is always minimised at s = q.\n')
else
    fprintf('\nmin U(q,s)/U(q,s-1) is achieved at s = \n')
    disp(argminratios)
end

fprintf(['\nMax. error in U ratios is no more than ' num2str(ratios_error) '.\n'])
C_error = C_vector.upper() - C_vector.lower();
disp(['Max. error in C is no more than ' num2str(max(C_error, [], 'all')) '.'])
    
C_output = C_vector.upper();

output1 = [qvalues', minratios', C_output', big_V_upper'];

output1(:, 2) = floor(output1(:, 2) * 1e6) / 1e6; % Makes ratios round down.
output1(:, [3 4]) = ceil(output1(:, [3 4]) * 1e6) / 1e6; % Makes everything else round up.
output1 = num2cell(output1.');
for j = 1 : numel(output1)
    if isnan(output1{j})
        output1{j} = [];
    end
end

fprintf('\nRelevant computed values follow:\n')

if output_latex
    disp('    q     Min U ratio            C(q)      V u. bound')
    fprintf('%5d & %13.6f & %13.6f & %13.6f \\\\\n', output1{:});
else
    disp('    q   Min U ratio          C(q)    V u. bound')
    fprintf('%5d %13.6f %13.6f %13.6f\n', output1{:});
end

ratios_increasing = all(minratios(1 : end - 1) < minratios(2 : end));
C_decreasing = all(C_output(1 : end - 1) > C_output(2 : end));
V_decreasing = all(V_upper(1 : end - 1) > V_upper(2 : end));

if ratios_increasing
    fprintf('\nMinimum U ratios are increasing.\n')
end
if C_decreasing
    disp('C is decreasing.')
end
if V_decreasing
    disp('Upper bound on V is decreasing.')
end

save("Covering_systems_data.mat")
dur = datetime - starttime;
fprintf('\nTotal computational time:\n')
disp(dur)
end