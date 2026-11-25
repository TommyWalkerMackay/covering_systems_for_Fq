function display_cov_sys (sys, latex)
% Usage: display_cov_sys (sys, latex)

% Inputs a (potential) covering system sys in cell format and displays a
% readable representation of sys. If latex is true, displays text to be
% copied to a LaTeX file.

if nargin < 1
    error('Not enough input arguments.')
end
if nargin < 2
    latex = false;
end

if ~iscell(sys)
    error('sys must be a cell array.')
end

nhyps = size(sys, 1);
ndims = size(sys, 2);

for j = 1 : numel(sys)
    sysj = sys{j};

    if ischar(sysj)
        if isscalar(sysj) && sysj == ':'
            sys{j} = '*';
        else
            error('Invalid entry in sys.')
        end
    elseif isnumeric(sysj)
        if isscalar(sysj)
            sys{j} = num2str(sysj);
        else
            error('Invalid entry in sys.')
        end
    else
        error('Invalid entry in sys.')
    end
end

d = num2cell((1 : ndims)');
n = num2cell((1 : nhyps)');
data = [n, sys];
data = data';

if latex

    fprintf(['\\#', repmat(' & $B_%d$', [1, ndims]), ' \\\\ \\hline\n'], d{:})
    fprintf(['%d', repmat(' & $%c$', [1, ndims]), ' \\\\\n'], data{:})

else

    fprintf([' #', repmat('   B_%d', [1, ndims]), '\n'], d{:})
    fprintf(['%2d ', repmat('   %c  ', [1, ndims]), '\n'], data{:})

end
end