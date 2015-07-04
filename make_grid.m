function make_grid()
%   Create a square grid of points mapping the interior of the domain
%   of the problem.  Each point contains 0 or 1.

%   First example: the unit circle.  Say that the grid is defined for all
%   points 0, +/-  0.1, +/- 0.2, ...

h = 0.1;   % grid spacing
N = 21;    % total number of rows/cols in our grid

grid = zeros(N,N);
count_interior_points = 0;

for i = 1:N
    for j = 1:N
        x = -1. + h * (i-1);
        y = -1. + h * (j-1);
        if x^2 + y^2 < 1
           grid(i,j) = 1; 
           count_interior_points = count_interior_points + 1;
        end
    end
end

% Create plot of points
% disp(grid);
disp(count_interior_points);

% vector containing x-coordin
X = zeros(count_interior_points, 1);
Y = zeros(count_interior_points, 1);
current_point = 1;
for i = 1:N
    for j = 1:N
        if grid(i,j) == 1
           X(current_point) = -1. + h * (i-1);
           Y(current_point) = -1. + h * (j-1);
           current_point = current_point + 1;
        end
    end
end
scatter(X,Y);

end

