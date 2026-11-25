function check_cov_systems (latex)
% Usage: check_cov_systems (latex)

% Checks that the covering systems given in Table 1 are valid.
% If latex is set to true, output is displayed so that it can be copied to
% a LaTeX file. latex = false by default.

if nargin < 1
    latex = false;
end

sys1 = { 1 , ':', ':';
        ':',  1 , ':';
        ':', ':',  1 ;
         2 ,  2 ,  2 
        };

sys2 = { 1 , ':', ':', ':';
        ':',  1 , ':', ':';
        ':', ':',  1 , ':';
        ':', ':', ':',  1 ;
         2 ,  2 , ':', ':';
         3 , ':',  2 , ':';
        ':',  3 ,  3 , ':';
         3 , ':', ':',  2 ;
        ':',  2 , ':',  3 ;
         2 ,  3 ,  2 , ':'
        };

switch 2
    case 1 % Use hand-constructed 29-hyperplane cover.
        sys3 = { 1 , ':', ':', ':', ':';
                ':',  2 , ':', ':', ':';
                ':', ':',  3 , ':', ':';
                ':', ':', ':',  4 , ':';
                ':', ':', ':', ':',  4 ;
                 2 ,  1 , ':', ':', ':';
                 3 , ':',  1 , ':', ':';
                ':',  3 ,  2 , ':', ':';
                 4 , ':', ':',  1 , ':';
                ':',  4 , ':',  2 , ':';
                ':', ':',  4 ,  3 , ':';
                 2 ,  3 ,  1 , ':', ':';
                 2 ,  4 , ':',  1 , ':';
                 3 , ':',  4 ,  1 , ':';
                ':',  3 ,  4 ,  2 , ':';
                 2 ,  3 ,  4 ,  1 , ':';
                 4 , ':', ':', ':',  1 ;
                ':',  1 , ':', ':',  2 ;
                ':', ':',  2 , ':',  3 ;
                ':', ':', ':',  3 ,  2 ;
                 4 ,  1 , ':', ':',  3 ;
                 3 , ':',  2 , ':',  1 ;
                ':',  3 ,  1 , ':',  3 ;
                 2 , ':', ':',  3 ,  1 ;
                ':',  4 , ':',  3 ,  3 ;
                ':', ':',  4 ,  2 ,  1 ;
                 3 ,  1 ,  4 , ':',  3 ;
                 4 ,  3 , ':',  2 ,  2 ;
                 3 , ':',  2 ,  1 ,  2 ;
                };
    case 2 % Use improved 26-hyperplane cover.
        sys3 = { 1 , ':', ':', ':', ':';
                ':',  1 , ':', ':', ':';
                ':', ':',  1 , ':', ':';
                ':', ':', ':',  1 , ':';
                ':', ':', ':', ':',  1 ;
                 2 ,  2 , ':', ':', ':';
                 2 , ':',  2 , ':', ':';
                ':',  3 ,  4 , ':', ':';
                 3 , ':', ':',  2 , ':';
                ':',  4 , ':',  3 , ':';
                ':', ':',  3 ,  4 , ':';
                 4 , ':', ':', ':',  4 ;
                ':',  2 , ':', ':',  3 ;
                ':', ':',  2 , ':',  2 ;
                ':', ':', ':',  2 ,  2 ;
                 2 ,  4 ,  4 , ':', ':';
                 3 ,  2 , ':',  3 , ':';
                 2 , ':',  3 ,  2 , ':';
                ':',  3 ,  3 ,  3 , ':';
                 4 ,  2 , ':', ':',  2 ;
                 3 , ':',  2 , ':',  4 ;
                ':',  3 ,  2 , ':',  3 ;
                 4 , ':', ':',  2 ,  3 ;
                ':',  4 , ':',  4 ,  3 ;
                ':', ':',  4 ,  4 ,  2 ;
                 3 , ':',  4 ,  4 ,  4 ;
                };
end

for q = 3 : 5
    if q == 3
        sys = sys1;
    elseif q == 4
        sys = sys2;
    else
        sys = sys3;
    end
    dims = (q - 1) * ones(1, q);

    if is_discrete_cov_sys(sys, dims)
        fprintf(['\nCovering system for [', num2str(q - 1), ']^', num2str(q), ':\n\n'])
        display_cov_sys(sys, latex)
    else
        disp(['Covering system for [', num2str(q - 1), ']^', num2str(q), ' is invalid.'])
        keyboard;
    end
end
end