classdef BigFloat
    % BIGFLOAT  Interval floating-point arithmetic with a huge exponent range.
    %
    % Each endpoint is represented as
    %
    %       m * 2^e
    %
    % where 0.5 <= abs(m) < 1, or m == 0.
    %
    % The object represents the interval [lo, hi].
    %
    % Arithmetic is outward rounded by moving the double mantissa one
    % conservative ULP after each elementary operation.
    %
    % This implementation is intended for finite values and ordinary
    % normal double mantissas.
    %
    % Much of this class was written using GPT 5.6 Luna.

    properties
        lo_m double
        lo_e int64

        hi_m double
        hi_e int64
    end

    methods

        %% Constructor
        function obj = BigFloat(x, e)
            if nargin == 0
                obj.lo_m = 0;
                obj.lo_e = int64(0);
                obj.hi_m = 0;
                obj.hi_e = int64(0);
                return
            end

            if nargin == 1
                if isa(x, 'BigFloat')
                    obj = x;
                    return
                end

                % Convert ordinary numeric input to an exact singleton
                % where possible.
                [m,k] = BigFloat.normalize(double(x));

                obj.lo_m = m;
                obj.lo_e = int64(k);
                obj.hi_m = m;
                obj.hi_e = int64(k);
                return
            end

            % Explicit mantissa/exponent constructor.
            [m,k] = BigFloat.normalize(double(x));

            obj.lo_m = m;
            obj.lo_e = int64(k) + int64(e);

            obj.hi_m = m;
            obj.hi_e = int64(k) + int64(e);
        end


        %% Addition
        function z = plus(a,b)
            [a,b] = BigFloat.promote(a,b);

            % Lower endpoint:
            % lower(a+b) = lower(a) + lower(b)
            [lm,le] = BigFloat.addLower( ...
                a.lo_m,a.lo_e,b.lo_m,b.lo_e);

            % Upper endpoint:
            % upper(a+b) = upper(a) + upper(b)
            [um,ue] = BigFloat.addUpper( ...
                a.hi_m,a.hi_e,b.hi_m,b.hi_e);

            % Ensures adding 0 has no effect.
            ais0 = BigFloat.isZero(a);
            bis0 = BigFloat.isZero(b);
            lm(ais0) = b.lo_m(ais0);
            le(ais0) = b.lo_e(ais0);
            um(ais0) = b.hi_m(ais0);
            ue(ais0) = b.hi_e(ais0);
            lm(bis0) = a.lo_m(bis0);
            le(bis0) = a.lo_e(bis0);
            um(bis0) = a.hi_m(bis0);
            ue(bis0) = a.hi_e(bis0);

            z = BigFloat.fromEndpoints(lm,le,um,ue);
        end


        %% Subtraction
        function z = minus(a,b)
            [a,b] = BigFloat.promote(a,b);

            % lower(a-b) = lower(a) - upper(b)
            [lm,le] = BigFloat.subLower( ...
                a.lo_m,a.lo_e,b.hi_m,b.hi_e);

            % upper(a-b) = upper(a) - lower(b)
            [um,ue] = BigFloat.subUpper( ...
                a.hi_m,a.hi_e,b.lo_m,b.lo_e);

            % Ensures subtracting 0 has no effect.
            ais0 = BigFloat.isZero(a);
            bis0 = BigFloat.isZero(b);
            lm(ais0) = -b.lo_m(ais0);
            le(ais0) = b.lo_e(ais0);
            um(ais0) = -b.hi_m(ais0);
            ue(ais0) = b.hi_e(ais0);
            lm(bis0) = a.lo_m(bis0);
            le(bis0) = a.lo_e(bis0);
            um(bis0) = a.hi_m(bis0);
            ue(bis0) = a.hi_e(bis0);

            z = BigFloat.fromEndpoints(lm,le,um,ue);
        end


        %% Unary minus
        function z = uminus(a)
            z = BigFloat.fromEndpoints( ...
                -a.hi_m,a.hi_e, ...
                -a.lo_m,a.lo_e);
        end


        %% Multiplication
        function z = times(a,b)
            [a,b] = BigFloat.promote(a,b);

            % All four endpoint products are considered.
            [m1,e1] = BigFloat.mulDirected( ...
                a.lo_m,a.lo_e,b.lo_m,b.lo_e,'down');

            [m2,e2] = BigFloat.mulDirected( ...
                a.lo_m,a.lo_e,b.hi_m,b.hi_e,'down');

            [m3,e3] = BigFloat.mulDirected( ...
                a.hi_m,a.hi_e,b.lo_m,b.lo_e,'down');

            [m4,e4] = BigFloat.mulDirected( ...
                a.hi_m,a.hi_e,b.hi_m,b.hi_e,'down');

            % For the lower endpoint we need the minimum.
            [lm,le] = BigFloat.minEndpoint( ...
                m1, e1, m2, e2, m3, e3, m4, e4);

            [u1,f1] = BigFloat.mulDirected( ...
                a.lo_m,a.lo_e,b.lo_m,b.lo_e,'up');

            [u2,f2] = BigFloat.mulDirected( ...
                a.lo_m,a.lo_e,b.hi_m,b.hi_e,'up');

            [u3,f3] = BigFloat.mulDirected( ...
                a.hi_m,a.hi_e,b.lo_m,b.lo_e,'up');

            [u4,f4] = BigFloat.mulDirected( ...
                a.hi_m,a.hi_e,b.hi_m,b.hi_e,'up');

            [um,ue] = BigFloat.maxEndpoint( ...
                u1, f1, u2, f2, u3, f3, u4, f4);

            % Ensures multiplying by 0 always gives 0.
            is0 = BigFloat.isZero(a) | BigFloat.isZero(b);
            lm(is0) = 0;
            le(is0) = 0;
            um(is0) = 0;
            ue(is0) = 0;

            z = BigFloat.fromEndpoints(lm,le,um,ue);
        end


        %% Division
        function z = rdivide(a,b)
            [a,b] = BigFloat.promote(a,b);

            if any(BigFloat.containsZeroInternal(b), 'all')
                error('BigFloat:DivisionByZero', ...
                    'Division interval contains zero.');
            end

            % Four endpoint quotients.
            [m1,e1] = BigFloat.divDirected( ...
                a.lo_m,a.lo_e,b.lo_m,b.lo_e,'down');

            [m2,e2] = BigFloat.divDirected( ...
                a.lo_m,a.lo_e,b.hi_m,b.hi_e,'down');

            [m3,e3] = BigFloat.divDirected( ...
                a.hi_m,a.hi_e,b.lo_m,b.lo_e,'down');

            [m4,e4] = BigFloat.divDirected( ...
                a.hi_m,a.hi_e,b.hi_m,b.hi_e,'down');

            [lm,le] = BigFloat.minEndpoint( ...
                m1, e1, m2, e2, m3, e3, m4, e4);

            [u1,f1] = BigFloat.divDirected( ...
                a.lo_m,a.lo_e,b.lo_m,b.lo_e,'up');

            [u2,f2] = BigFloat.divDirected( ...
                a.lo_m,a.lo_e,b.hi_m,b.hi_e,'up');

            [u3,f3] = BigFloat.divDirected( ...
                a.hi_m,a.hi_e,b.lo_m,b.lo_e,'up');

            [u4,f4] = BigFloat.divDirected( ...
                a.hi_m,a.hi_e,b.hi_m,b.hi_e,'up');

            [um,ue] = BigFloat.maxEndpoint( ...
                u1, f1, u2, f2, u3, f3, u4, f4);

            z = BigFloat.fromEndpoints(lm,le,um,ue);
        end


        %% Multiplication operator
        function z = mtimes(a,b) 
            [a,b] = BigFloat.promote(a,b); 
            if ~ismatrix(a) || ~ismatrix(b) 
                error('Arguments must be 2-dimensional, or at least one argument must be scalar.') 
            end 
            if isscalar(a) || isscalar(b)
                z = a .* b; 
                return 
            end 
            if size(a, 2) ~= size(b, 1) 
                error('Incompatible dimensions for matrix multiplication.') 
            end 

            alm = a.lo_m;
            ale = a.lo_e;
            ahm = a.hi_m;
            ahe = a.hi_e;
            blm = b.lo_m;
            ble = b.lo_e;
            bhm = b.hi_m;
            bhe = b.hi_e;

            szz = [size(a, 1), size(b, 2)];

            zlm = zeros(szz);
            zle = int64(zeros(szz));
            zhm = zeros(szz);
            zhe = int64(zeros(szz));
 
            for j = 1 : size(a, 2)
                almj = alm(:, j);
                alej = ale(:, j);
                ahmj = ahm(:, j);
                ahej = ahe(:, j);
                blmj = blm(j, :);
                blej = ble(j, :);
                bhmj = bhm(j, :);
                bhej = bhe(j, :);

                % All four endpoint products are considered.
                [m1,e1] = BigFloat.mulDirected( ...
                    almj,alej,blmj,blej,'down');

                [m2,e2] = BigFloat.mulDirected( ...
                    almj,alej,bhmj,bhej,'down');

                [m3,e3] = BigFloat.mulDirected( ...
                    ahmj,ahej,blmj,blej,'down');

                [m4,e4] = BigFloat.mulDirected( ...
                    ahmj,ahej,bhmj,bhej,'down');

                % For the lower endpoint we need the minimum.
                [lm,le] = BigFloat.minEndpoint( ...
                    m1, e1, m2, e2, m3, e3, m4, e4);

                [u1,f1] = BigFloat.mulDirected( ...
                    almj,alej,blmj,blej,'up');

                [u2,f2] = BigFloat.mulDirected( ...
                    almj,alej,bhmj,bhej,'up');

                [u3,f3] = BigFloat.mulDirected( ...
                    ahmj,ahej,blmj,blej,'up');

                [u4,f4] = BigFloat.mulDirected( ...
                    ahmj,ahej,bhmj,bhej,'up');

                [um,ue] = BigFloat.maxEndpoint( ...
                    u1, f1, u2, f2, u3, f3, u4, f4);

                % Ensures multiplying by 0 always gives 0.
                is0 = (almj == 0 & ahmj == 0) | (blmj == 0 & bhmj == 0);
                lm(is0) = 0;
                le(is0) = 0;
                um(is0) = 0;
                ue(is0) = 0;

                [zlm, zle] = BigFloat.addLower(zlm, zle, lm, le);
                [zhm, zhe] = BigFloat.addUpper(zhm, zhe, um, ue);
                
            end 
            z = BigFloat.fromEndpoints(zlm, zle, zhm, zhe);
        end


        %% Division operator
        function z = mrdivide(a,b)
            if isscalar(b)
                z = a ./ b;
            else
                error('Matrix division not yet implemented.')
            end
        end


        %% Sum
        function y = sum(x, dim)

            if nargin < 2
                dim = find(size(x.lo_m) ~= 1, 1);
                if isempty(dim)
                    dim = 1;
                end
            end

            [lm, le] = BigFloat.sumLower(x.lo_m, x.lo_e, dim);
            [um, ue] = BigFloat.sumUpper(x.hi_m, x.hi_e, dim);

            y = BigFloat.fromEndpoints(lm, le, um, ue);
        end


        %% Prod
        function y = prod(x, dim)

            if nargin < 2
                dim = find(size(x.lo_m) ~= 1, 1);
                if isempty(dim)
                    dim = 1;
                end
            end

            [lm, le] = BigFloat.prodLower(x.lo_m, x.lo_e, dim);
            [um, ue] = BigFloat.prodUpper(x.hi_m, x.hi_e, dim);

            y = BigFloat.fromEndpoints(lm, le, um, ue);
        end


        %% Cumulative sum
        function y = cumsum(x, dim)

            if nargin < 2
                dim = find(size(x.lo_m) ~= 1, 1);
                if isempty(dim)
                    dim = 1;
                end
            end

            [lm, le] = BigFloat.cumsumLower(x.lo_m, x.lo_e, dim);
            [um, ue] = BigFloat.cumsumUpper(x.hi_m, x.hi_e, dim);

            y = BigFloat.fromEndpoints(lm, le, um, ue);
        end


        %% Cumulative product
        function y = cumprod(x, dim)

            if nargin < 2
                dim = find(size(x.lo_m) ~= 1, 1);
                if isempty(dim)
                    dim = 1;
                end
            end

            [lm, le] = BigFloat.cumprodLower(x.lo_m, x.lo_e, dim);
            [um, ue] = BigFloat.cumprodUpper(x.hi_m, x.hi_e, dim);

            y = BigFloat.fromEndpoints(lm, le, um, ue);
        end

        %% Combine the range of two BigFloats
        function z = combine_range(a,b)
            [a,b] = BigFloat.promote(a,b);
            [lm, le] = BigFloat.minEndpoint(a.lo_m, a.lo_e, b.lo_m, b.lo_e);
            [um, ue] = BigFloat.maxEndpoint(a.hi_m, a.hi_e, b.hi_m, b.hi_e);
            z = BigFloat.fromEndpoints(lm, le, um, ue);
        end


        %% Display 
        function disp(x)
            fprintf('Lower bound:\n'); 
            fprintf(' mantissa =\n'); 
            disp(x.lo_m); 
            fprintf(' exponent =\n'); 
            disp(x.lo_e); 
            fprintf('Upper bound:\n'); 
            fprintf(' mantissa =\n'); 
            disp(x.hi_m); 
            fprintf(' exponent =\n'); 
            disp(x.hi_e); 
        end


        %% Convert to double
        function y = double(x)
            lo = x.lo_m .* 2.^double(x.lo_e);
            hi = x.hi_m .* 2.^double(x.hi_e);
            if all(lo == hi, 'all') 
                y = lo; 
            else 
                y = [lo(:); hi(:)]; 
                y = reshape(y, [size(lo), 2]); 
            end 
        end


        %% Return lower endpoint as a double where representable
        function y = lower(x)
            y = x.lo_m .* 2.^double(x.lo_e);
        end


        %% Return upper endpoint as a double where representable
        function y = upper(x)
            y = x.hi_m .* 2.^double(x.hi_e);
        end


        %% Normal property assignment
        function obj = subsasgn(obj, S, value)

            if strcmp(S(1).type, '()')
        
                idx = S(1).subs;
        
                % Convert RHS to BigFloat if necessary
                if ~isa(value, 'BigFloat')
                    value = BigFloat(value);
                end
        
                obj.lo_m(idx{:}) = value.lo_m;
                obj.lo_e(idx{:}) = value.lo_e;
                obj.hi_m(idx{:}) = value.hi_m;
                obj.hi_e(idx{:}) = value.hi_e;
        
                % Handle further indexing, e.g. S(1,2,3).lo_m
                if numel(S) > 1
                    obj = builtin('subsasgn', obj, S(2:end), value);
                end
        
            elseif strcmp(S(1).type, '.')
        
                obj = builtin('subsasgn', obj, S, value);
        
            else
                error('Unsupported subscript type.');
            end
        
        end


        %% Normal element access
        function out = subsref(obj, S)

            out = BigFloat.subsrefInternal(obj, S);

        end


        %% Return size
        function varargout = size(obj, dim)

            if nargin == 1
                [varargout{1:nargout}] = size(obj.lo_m);
            else
                varargout{1} = size(obj.lo_m, dim);
            end

        end


        %% Return numel
        function n = numel(obj)
            n = numel(obj.lo_m);
        end


        %% Check if the object is empty
        function tf = isempty(obj)
            tf = isempty(obj.lo_m);
        end

        %% Transpose
        function obj = transpose(obj)
            obj.lo_m = obj.lo_m.';
            obj.lo_e = obj.lo_e.';
            obj.hi_m = obj.hi_m.';
            obj.hi_e = obj.hi_e.';
        end

        %% Transpose
        function obj = ctranspose(obj)
            obj.lo_m = obj.lo_m.';
            obj.lo_e = obj.lo_e.';
            obj.hi_m = obj.hi_m.';
            obj.hi_e = obj.hi_e.';
        end


        %% Reshape
        function obj = reshape(obj, varargin)
            obj.lo_m = reshape(obj.lo_m, varargin{:});
            obj.lo_e = reshape(obj.lo_e, varargin{:});
            obj.hi_m = reshape(obj.hi_m, varargin{:});
            obj.hi_e = reshape(obj.hi_e, varargin{:});
        end


        %% Concatenation
        function obj = cat(dim, varargin)

            % Promote numeric inputs
            for k = 1:numel(varargin)
                if ~isa(varargin{k}, 'BigFloat')
                    varargin{k} = BigFloat(varargin{k});
                end
            end
    
            lm = cellfun(@(x) x.lo_m, varargin, 'UniformOutput', false);
            le = cellfun(@(x) x.lo_e, varargin, 'UniformOutput', false);
            um = cellfun(@(x) x.hi_m, varargin, 'UniformOutput', false);
            ue = cellfun(@(x) x.hi_e, varargin, 'UniformOutput', false);
    
            obj = BigFloat.fromEndpoints( ...
                builtin('cat', dim, lm{:}), ...
                builtin('cat', dim, le{:}), ...
                builtin('cat', dim, um{:}), ...
                builtin('cat', dim, ue{:}));
        end


        %% Horizontal concatenation
        function obj = horzcat(varargin)
            obj = BigFloat.catInternal(2, varargin{:});
        end


        %% Vertical concatenation
        function obj = vertcat(varargin)
            obj = BigFloat.catInternal(1, varargin{:});
        end


        %% Length
        function n = length(obj)
            n = length(obj.lo_m);
        end


        %% End
        function ind = end(obj, k, n)

            if n == 1
                ind = numel(obj.lo_m);
            else
                sz = size(obj.lo_m);
        
                if k <= numel(sz)
                    ind = sz(k);
                else
                    ind = 1;
                end
            end
        
        end


        %% Check whether interval contains zero
        function tf = containsZero(x)

            tf = BigFloat.containsZeroInternal(x);
        end

    end


    methods (Static, Access = private)

        %% Promote ordinary numbers
        function [a,b] = promote(a,b)
            if ~isa(a,'BigFloat')
                a = BigFloat(a);
            end
            if ~isa(b,'BigFloat')
                b = BigFloat(b);
            end
        end


        %% Create directly from normalized endpoints
        function z = fromEndpoints(lm,le,um,ue)
            z = BigFloat();
            z.lo_m = lm;
            z.lo_e = int64(le);
            z.hi_m = um;
            z.hi_e = int64(ue);
        end        


        %% Normalize a double mantissa
        function [m,e] = normalize(x)

            [m,e] = log2(x);
        end


        %% Outward rounding of a normalized mantissa
        function [m,e] = roundMantissa(m,e,direction)

            % Normalize first.
            [m0,k] = log2(m);
            e = int64(e) + int64(k);
            m = m0;

            if strcmp(direction,'down')
                m = BigFloat.prevDouble(m);
            else
                m = BigFloat.nextDouble(m);
            end

            % Renormalize after the one-ULP movement.
            [m,k] = log2(m);
            e = e + int64(k);
        end


        %% Previous representable double.
        %
        % MATLAB does not provide nextafter, so construct the predecessor
        % using the IEEE-754 bit representation.
        function y = prevDouble(x)

            sz = size(x);
            x = x(:);
            xzero = x == 0;
            xpos = x > 0;
            xneg = x < 0;
            
            u = typecast(x,'uint64');

            u(xpos) = u(xpos) - uint64(1);
            u(xneg) = u(xneg) + uint64(1);
            
            y = typecast(u,'double');
            y(xzero) = -realmin;
            y = reshape(y, sz);
        end


        %% Next representable double.
        function y = nextDouble(x)

            sz = size(x);
            x = x(:);
            xzero = x == 0;
            xpos = x > 0;
            xneg = x < 0;

            u = typecast(x,'uint64');

            u(xpos) = u(xpos) + uint64(1);
            u(xneg) = u(xneg) - uint64(1);
            
            y = typecast(u,'double');
            y(xzero) = realmin;
            y = reshape(y, sz);
        end


        %% Addition rounded downward 
        function [m,e] = addLower(a,ae,b,be) 
            azero = a == 0; 
            bzero = b == 0; 
            E = max(ae,be); 
            da = double(ae-E); 
            db = double(be-E); 
            aa = a .* 2.^da; 
            bb = b .* 2.^db; 
            r = aa + bb; 
            [m,e] = BigFloat.roundMantissa(r,E,'down'); 
            [m(azero),e(azero)] = BigFloat.roundMantissa(b(azero),be(azero),'down'); 
            [m(bzero),e(bzero)] = BigFloat.roundMantissa(a(bzero),ae(bzero),'down'); 
        end 

        
        %% Addition rounded upward 
        function [m,e] = addUpper(a,ae,b,be) 
            azero = a == 0; 
            bzero = b == 0; 
            E = max(ae,be); 
            da = double(ae-E); 
            db = double(be-E); 
            aa = a .* 2.^da; 
            bb = b .* 2.^db; 
            r = aa + bb; 
            [m,e] = BigFloat.roundMantissa(r,E,'up'); 
            [m(azero),e(azero)] = BigFloat.roundMantissa(b(azero),be(azero),'up'); 
            [m(bzero),e(bzero)] = BigFloat.roundMantissa(a(bzero),ae(bzero),'up'); 
        end


        %% Subtraction rounded downward
        function [m,e] = subLower(a,ae,b,be)

            E = max(ae,be);

            da = double(ae-E);
            db = double(be-E);

            aa = a .* 2.^da;
            bb = b .* 2.^db;

            r = aa - bb;

            [m,e] = BigFloat.roundMantissa(r,E,'down');
        end


        %% Subtraction rounded upward
        function [m,e] = subUpper(a,ae,b,be)

            E = max(ae,be);

            da = double(ae-E);
            db = double(be-E);

            aa = a .* 2.^da;
            bb = b .* 2.^db;

            r = aa - bb;

            [m,e] = BigFloat.roundMantissa(r,E,'up');
        end


        %% Directed multiplication
        function [m,e] = mulDirected(a,ae,b,be,direction)

            r = a .* b;
            E = ae + be;

            [m,e] = BigFloat.roundMantissa(r,E,direction);
        end


        %% Directed division
        function [m,e] = divDirected(a,ae,b,be,direction)

            r = a ./ b;
            E = ae - be;

            [m,e] = BigFloat.roundMantissa(r,E,direction);
        end


        %% Lower sum
        function [m,e] = sumLower(m,e,dim)

            n = size(m,dim);

            idx = repmat({':'}, 1, ndims(m));
            idx{dim} = 1;

            sm = m(idx{:});
            se = e(idx{:});

            for k = 2:n
                idx{dim} = k;

                [sm,se] = BigFloat.addLower( ...
                    sm, se, m(idx{:}), e(idx{:}));
            end

            m = sm;
            e = se;
        end


        %% Upper sum
        function [m,e] = sumUpper(m,e,dim)

            n = size(m,dim);

            idx = repmat({':'}, 1, ndims(m));
            idx{dim} = 1;

            sm = m(idx{:});
            se = e(idx{:});

            for k = 2:n
                idx{dim} = k;

                [sm,se] = BigFloat.addUpper( ...
                    sm, se, m(idx{:}), e(idx{:}));
            end

            m = sm;
            e = se;
        end

           
        %% Lower cumulative sum
        function [m,e] = cumsumLower(m,e,dim)

            n = size(m,dim);

            for k = 2:n
                idx1 = repmat({':'},1,ndims(m));
                idx2 = idx1;

                idx1{dim} = k;
                idx2{dim} = k-1;

                [m(idx1{:}), e(idx1{:})] = ...
                    BigFloat.addLower( ...
                    m(idx1{:}), e(idx1{:}), ...
                    m(idx2{:}), e(idx2{:}));
            end
        end


        %% Upper cumulative sum
        function [m,e] = cumsumUpper(m,e,dim)

            n = size(m,dim);

            for k = 2:n
                idx1 = repmat({':'},1,ndims(m));
                idx2 = idx1;

                idx1{dim} = k;
                idx2{dim} = k-1;

                [m(idx1{:}), e(idx1{:})] = ...
                    BigFloat.addUpper( ...
                    m(idx1{:}), e(idx1{:}), ...
                    m(idx2{:}), e(idx2{:}));
            end
        end


        %% Lower cumulative product
        function [m,e] = cumprodLower(m,e,dim)

            n = size(m,dim);

            for k = 2:n
                idx1 = repmat({':'},1,ndims(m));
                idx2 = idx1;

                idx1{dim} = k;
                idx2{dim} = k-1;

                [m(idx1{:}), e(idx1{:})] = ...
                    BigFloat.mulDirected( ...
                    m(idx1{:}), e(idx1{:}), ...
                    m(idx2{:}), e(idx2{:}), 'down');
            end
        end


        %% Upper cumulative product
        function [m,e] = cumprodUpper(m,e,dim)

            n = size(m,dim);

            for k = 2:n
                idx1 = repmat({':'},1,ndims(m));
                idx2 = idx1;

                idx1{dim} = k;
                idx2{dim} = k-1;

                [m(idx1{:}), e(idx1{:})] = ...
                    BigFloat.mulDirected( ...
                    m(idx1{:}), e(idx1{:}), ...
                    m(idx2{:}), e(idx2{:}), 'up');
            end
        end

        %% Check whether interval contains zero
        function tf = containsZeroInternal(x)

            tf = (x.lo_m == 0) | ...
                 (x.hi_m == 0) | ...
                 (x.lo_m < 0 & x.hi_m > 0);
        end


        %% Check which elements of x are zero.
        function tf = isZero(x)
            tf = x.lo_m == 0 & x.hi_m == 0;
        end

        %% Find minimum endpoint
        function [m,e] = minEndpoint(m1,e1,m2,e2,m3,e3,m4,e4)

            if nargin < 2 || mod(nargin, 2) ~= 0 || nargin > 8
                error('Incorrect number of inputs');
            end
            k_max = nargin / 2;

            mbest = m1;
            ebest = e1;

            for k = 2:k_max
                if k == 2
                    mthis = m2;
                    ethis = e2;
                elseif k == 3
                    mthis = m3;
                    ethis = e3;
                else
                    mthis = m4;
                    ethis = e4;
                end
                better = BigFloat.lessThan(mthis,ethis,mbest,ebest);
                mbest(better) = mthis(better);
                ebest(better) = ethis(better);
                
            end

            m = mbest;
            e = ebest;
        end


        %% Find maximum endpoint
        function [m,e] = maxEndpoint(m1,e1,m2,e2,m3,e3,m4,e4)

            if nargin < 2 || mod(nargin, 2) ~= 0 || nargin > 8
                error('Incorrect number of inputs');
            end
            k_max = nargin / 2;

            mbest = m1;
            ebest = e1;

            for k = 2:k_max
                if k == 2
                    mthis = m2;
                    ethis = e2;
                elseif k == 3
                    mthis = m3;
                    ethis = e3;
                else
                    mthis = m4;
                    ethis = e4;
                end
                better = BigFloat.lessThan(mbest,ebest,mthis,ethis);
                mbest(better) = mthis(better);
                ebest(better) = ethis(better);
                
            end

            m = mbest;
            e = ebest;
        end


        %% Compare two normalized floating-point values
        function tf = lessThan(a,ae,b,be)

            sgn_less = sign(a) < sign(b);
            sgn_pos = sign(a) == 1 & sign(b) == 1;
            sgn_neg = sign(a) == -1 & sign(b) == -1;

            tf = sgn_less | (sgn_pos & (ae < be | (ae == be & a < b)))...
                | sgn_neg & (ae > be | (ae == be & a < b));

        end


        %% Concatenation
        function obj = catInternal(dim, varargin)

            % Promote numeric arguments
            for k = 1:numel(varargin)
                if ~isa(varargin{k}, 'BigFloat')
                    varargin{k} = BigFloat(varargin{k});
                end
            end

            % Extract underlying arrays
            lm = cell(1,numel(varargin));
            le = cell(1,numel(varargin));
            hm = cell(1,numel(varargin));
            he = cell(1,numel(varargin));

            for k = 1:numel(varargin)
                lm{k} = varargin{k}.lo_m;
                le{k} = varargin{k}.lo_e;
                hm{k} = varargin{k}.hi_m;
                he{k} = varargin{k}.hi_e;
            end

            % Concatenate the arrays
            obj = BigFloat.fromEndpoints( ...
                builtin('cat',dim,lm{:}), ...
                builtin('cat',dim,le{:}), ...
                builtin('cat',dim,hm{:}), ...
                builtin('cat',dim,he{:}));
        end


        %% Normal element access
        function out = subsrefInternal(obj, S)

            if strcmp(S(1).type, '()')

                idx = S(1).subs;

                out = BigFloat();

                out.lo_m = obj.lo_m(idx{:});
                out.lo_e = obj.lo_e(idx{:});
                out.hi_m = obj.hi_m(idx{:});
                out.hi_e = obj.hi_e(idx{:});

                % MATLAB removes trailing singleton dimensions for numeric arrays.
                % Force the same behaviour for the BigFloat wrapper.
                sz = size(out.lo_m);
                if numel(sz) > 2
                    while numel(sz) > 2 && sz(end) == 1
                        sz(end) = [];
                    end

                    out.lo_m = reshape(out.lo_m, sz);
                    out.lo_e = reshape(out.lo_e, sz);
                    out.hi_m = reshape(out.hi_m, sz);
                    out.hi_e = reshape(out.hi_e, sz);
                end

                if numel(S) > 1
                    out = subsref(out, S(2:end));
                end

            elseif strcmp(S(1).type, '.')

                out = builtin('subsref', obj, S);

            else
                error('Unsupported subscript type.');
            end

        end
    end
end