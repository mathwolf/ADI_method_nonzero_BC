function stencil = make_stencil(h, N, x_min, y_min)
%   Create a square grid of points mapping the interior of the domain
%   of the problem.  Each point contains 0 or 1.
%   h: grid spacing
%   N: total number of rows and cols in the grid
%   x_min: starting value for x
%   y_min: ... for y

%   First example: the unit circle.  Say that the grid is defined for all
%   points 0, +/-  0.1, +/- 0.2, ...

stencil = zeros(N,N);
% n_pts = 0;

for i = 1:N
    for j = 1:N
        x = x_min + h * (i-1);
        y = y_min + h * (j-1);
        if (phi1(x) < y) && (y < phi2(x)) && ...
                (psi1(y) < x) && (x < psi2(y))
           stencil(i,j) = 1; 
           % n_pts = n_pts + 1;
        end
    end
end

%{
% Create plot of points
disp(n_pts);

% vector containing x-coordinate
X = zeros(n_pts, 1);
Y = zeros(n_pts, 1);
current_point = 1;
for i = 1:N
    for j = 1:N
        if stencil(i,j) == 1
           X(current_point) = -1. + h * (i-1);
           Y(current_point) = -1. + h * (j-1);
           current_point = current_point + 1;
        end
    end
end
scatter(X,Y);
%}
end
