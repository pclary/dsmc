function [xn, yn] = rk4step(t, h, x, y, u, v, varargin)
%RK4STEP Steps 2D points forward in time according to the velocity field
%specified by (u, v) using a 4th order Runge-Kutta method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = u(t, x, y, varargin{:});
v1 = v(t, x, y, varargin{:});
x1 = x + h/2 * u1;
y1 = y + h/2 * v1;
u2 = u(t + h/2, x1, y1, varargin{:});
v2 = v(t + h/2, x1, y1, varargin{:});
x2 = x + h/2 * u2;
y2 = y + h/2 * v2;
u3 = u(t + h/2, x2, y2, varargin{:});
v3 = v(t + h/2, x2, y2, varargin{:});
x3 = x + h * u3;
y3 = y + h * v3;
u4 = u(t + h, x3, y3, varargin{:});
v4 = v(t + h, x3, y3, varargin{:});
xn = x + h/6 * (u1 + 2 * u2 + 2 * u3 + u4);
yn = y + h/6 * (v1 + 2 * v2 + 2 * v3 + v4);

end

