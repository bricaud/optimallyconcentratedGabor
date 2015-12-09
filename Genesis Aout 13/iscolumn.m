% Copyright (C) 1996, 1997, 2002, 2004, 2005, 2006, 2007, 2008 John W. Eaton
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {Function File} {} iscolumn (@var{a})
% Return 1 if @var{a} is a column vector.  Otherwise, return 0.
% @seealso{size, rows, columns, length, isscalar, ismatrix, isvector, isrow}
% @end deftypefn

% Author: Carlos Deambroggio (2011)
% Author of isvector(): jwe

function retval = iscolumn (x)

  retval = 0;

  if (nargin == 1)
    sz = size (x);
    retval = (ndims (x) == 2 && (sz(2) == 1));
  else
    print_usage ();
  end

end

%!assert(iscolumn (1));

%!assert(iscolumn ([1; 2; 3]));

%!assert(!(iscolumn ([])));

%!assert((iscolumn ([12; 34])));

%!assert(!(iscolumn ([13 35])));

%!test
%! warning("off","Octave:single-quote-string");
%! assert((iscolumn ("t")));

%!test
%! warning("off","Octave:single-quote-string");
%! assert((iscolumn (["t";"e";"s";"t"])));

%!test
%! warning("off","Octave:single-quote-string");
%! assert((iscolumn ([["t";"e";"s";"t"]; ["i";"n";"g"]])));

%!assert(iscolumn(zeros(0,1)));

%!assert(!iscolumn(zeros(0,18)));

%!test
%! s.a = 1;
%! assert((iscolumn (s)));

%!error iscolumn ();

%!error iscolumn ([1, 2], 2);

