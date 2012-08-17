function ge = gb2ge (varargin)
% GB2GE Converts LAPACK-like general banded matrices to dense storage.
% The general band matrix gb must be stored following LAPACK conventions.
% m is the number of rows, n the number of columns, and kl and ku
% are the number of sub- and super-diagonals, respectively.
%
% Invoke as any of
%          ge = gb2ge (gb, m, n, kl, ku, ld)
%          ge = gb2ge (gb, m, n, kl, ku)
%          ge = gb2ge (gb, n, kl, ku)
%          ge = gb2ge (gb, n, k)
% where the first form is for general rectangular matrices, the second
% for general square matrices, and the last for square matrices with
% the same number of sub- and super-diagonals.

switch (nargin)
    case 6
        gb = varargin{1};
        m  = varargin{2};
        n  = varargin{3};
        kl = varargin{4};
        ku = varargin{5};
        ld = varargin{6};
    case 5
        gb = varargin{1};
        m  = varargin{2};
        n  = varargin{3};
        kl = varargin{4};
        ku = varargin{5};
        ld = kl + 1 + ku;
    case 4
        gb = varargin{1};
        m  = varargin{2};
        n  = varargin{2};
        kl = varargin{3};
        ku = varargin{4};
        ld = kl + 1 + ku;
    case 3
        gb = varargin{1};
        m  = varargin{2};
        n  = varargin{2};
        kl = varargin{3};
        ku = varargin{3};
        ld = kl + 1 + ku;
    otherwise
        error(nargchk(3, 6, nargin));
end

ge = zeros(m, n);
for j = 1 : m
    for i = max(1, j - ku) : min(m, j + kl)
        ge(i, j) = gb((j-1)*ld + ku + i - j + 1);
    end
end

end

%!test
%! % From http://www.netlib.org/lapack/lug/node124.html
%! m  = 5;
%! n  = 5;
%! kl = 2;
%! ku = 1;
%! gb = [ NaN, 11,  21,  31, ...
%!         12, 22,  32,  42, ...
%!         23, 33,  43,  53, ...
%!         34, 44,  54, NaN, ...
%!         45, 55, NaN, NaN];
%! ge = [ 11, 12,  0,  0,  0 ; ...
%!        21, 22, 23,  0,  0 ; ...
%!        31, 32, 33, 34,  0 ; ...
%!         0, 42, 43, 44, 45 ; ...
%!         0,  0, 53, 54, 55 ];
%! assert(all(ge == gb2ge(gb, m, n, kl, ku, kl + 1 + ku)));
%! assert(all(ge == gb2ge(gb, m, n, kl, ku)));
%! assert(all(ge == gb2ge(gb, n, kl, ku)));
