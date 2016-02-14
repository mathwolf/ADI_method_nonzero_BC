function tester()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

TRUE = 1;
FALSE = 0;

% Constants used to switch between different test functions.
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLYNOMIAL = 4;
global test_solution
test_solution = EXPONENT_1;

% Constants used to switch between different test domains.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
global domain
domain = CIRCLE;

% Description of spatial grid.  We divide both the x and y dimensions of
% the problem into the same number of gridpoints.  First, identify the min 
% and max values in both directions.  The actual domain of the PDE will be 
% a subset of the big rectangle defined by these values.

% For now, we are just defaulting to a cube of side length 2.
if domain == ELLIPSE
   x_min = -1;
   x_max = 1;
   y_min = -0.5;
   y_max = 0.5;
elseif domain == RECTANGLE
   x_min = -1;
   x_max = 1;
   y_min = -0.5;
   y_max = 0.5;    
elseif domain == DIAMOND_2
   x_min = -1;
   x_max = 0.5;
   y_min = -1;
   y_max = 0.5;    
else
    % all the other test domains have the same spatial range
   x_min = -1;
   x_max = 1;
   y_min = -1;
   y_max = 1;
   z_min = -1;
   z_max = 1;
end

N = 9;
hx = (x_max - x_min)/N;
hy = (y_max - y_min)/N;
hz = (z_max - z_min)/N;

% Define the gridpoints that are part of the interior of the domain of the
% PDE.  Check each gridpoint in the array to see if it is an interior 
% point.  Examine one stack at a time beginning with the leftmost column
% and uppermost row.
for i = 1:N-1
   for j = 1:N-1
       for k = 1:N-1
           x = x_min + hx * i;
           y = y_min + hy * j;
           z = z_min + hz * k;
           if (phi1(x) < y) && (y < phi2(x)) && ...
                    (zeta1(x,y) < z) && (z < zeta2(x,y))
               % The point is an interior point
               grid(i,j,k).on = TRUE; 
               grid(i,j,k).x = x;
               grid(i,j,k).y = y;
               grid(i,j,k).z = z;
               grid(i,j,k).U = g3(x,y,z);   % use initial condition to fill in the
                                            % value of the approximation
               grid(i,j,k).V = 0.;          % grid function for partial timesteps
           else
               % The point is an exterior point
               grid(i,j,k).on = FALSE;
           end
       end
   end
end

table = zeros(N-1,N-1,N-1);
% Print for debugging.
for i = 1:N-1
    for j = 1:N-1
        for k = 1:N-1
            table(i,j,k) = grid(i,j,k).on;
        end
    end
end

disp(table);

% Create a data structure that describes the rows of the grid.
n_rows = 0;
for j = 1:N-1
    for k = 1:N-1
        % Go along the row until we reach the end or find an interior point.
        i = 1;
        while (i <= N-1) && (grid(i,j,k).on == FALSE)
           i = i + 1; 
        end

        if i > N-1
           % We are at the end of an empty row, go to the next row
           continue;
        end

        % We are at the first interior point in a new row
        n_rows = n_rows + 1;    
        row(n_rows).j = j;      % assign the next number in sequence as the
        row(n_rows).k = k;      % index for the current row
        row(n_rows).i_min = i;

        % Go along the row until we reach the end or or find an exterior point
        while (i <= N-1) && (grid(i,j,k).on == TRUE)
           i = i + 1; 
        end
        row(n_rows).i_max = i - 1;

        % Add boundary terms for first and last point in current row
        row(n_rows).btf.y = y_min + hy * j;
        row(n_rows).btf.z = z_min + hz * k;
        row(n_rows).btf.x = eta1(row(n_rows).btf.y, row(n_rows).btf.z);
        row(n_rows).btf.h_prime = x_min + hx * row(n_rows).i_min ...
            - row(n_rows).btf.x;
        row(n_rows).btf.U = g3(row(n_rows).btf.x, row(n_rows).btf.y, ...
            row(n_rows).btf.z);
                % using initial conditions for first value of U

        row(n_rows).btl.y = y_min + hy * j;
        row(n_rows).btl.z = z_min + hz * k;
        
        row(n_rows).btl.x = eta2(row(n_rows).btl.y, row(n_rows).btl.z);
        row(n_rows).btl.h_prime = row(n_rows).btl.x ...
            -  ( x_min + hx * row(n_rows).i_max );   
        row(n_rows).btl.U = g3(row(n_rows).btl.x, row(n_rows).btl.y, ...
            row(n_rows).btl.z);
                % using initial conditions
    end
end

% Create a data structure that describes the columns of the grid.
n_cols = 0;
for i = 1:N-1
    for k = 1:N-1
        % Go through col until we reach the end or find an interior point
        j = 1;
        while (j <= N - 1) && (grid(i,j,k).on == FALSE) 
           j = j + 1; 
        end

        if j > N - 1
           % We are at the end of an empty col, go to the next col
           continue;
        end

        % We are at the first interior point in a new col
        n_cols = n_cols + 1;
        col(n_cols).i = i;      % assign the next number in sequence as the
        col(n_cols).k = k;      % index for the current column
        col(n_cols).j_min = j;

        % Go through the column until we reach the end or find an exterior
        % point
        while (j <= N-1) && (grid(i,j,k).on == TRUE) 
           j = j + 1; 
        end
        col(n_cols).j_max = j - 1;

        % Add boundary terms for first and last points in current column
        col(n_cols).btf.x = x_min + hx * i;
        col(n_cols).btf.z = z_min + hz * k;
        col(n_cols).btf.y = theta1(col(n_cols).btf.x, col(n_cols).btf.z);
        col(n_cols).btf.h_prime = y_min + hy * col(n_cols).j_min ...
            - col(n_cols).btf.y;    
        col(n_cols).btf.U = g3(col(n_cols).btf.x, col(n_cols).btf.y, ...
            col(n_cols).btf.z);
                % use initial conditions for first value of U

        col(n_cols).btl.x = x_min + hx * i;
        col(n_cols).btl.z = z_min + hz * k;
        col(n_cols).btl.y = theta2(col(n_cols).btl.x, col(n_cols).btl.z);
        col(n_cols).btl.h_prime = col(n_cols).btl.y ...
            - ( y_min + hy * col(n_cols).j_max  );
        col(n_cols).btl.U = g3(col(n_cols).btl.x, col(n_cols).btl.y, ...
            col(n_cols).btl.z);
                % use initial conditions
    end
end

% Create a data structure that describes the stacks of the grid.
n_stacks = 0;
for i = 1:N-1
    for j = 1:N-1
        % Go through stack until we reach the end or find an interior point
        k = 1;
        while (k <= N - 1) && (grid(i,j,k).on == FALSE) 
           k = k + 1; 
        end

        if k > N - 1
           % We are at the end of an empty stack, go to the next stack
           continue;
        end

        % We are at the first interior point in a new stack
        n_stacks = n_stacks + 1;
        stack(n_stacks).i = i;      % assign the next number in sequence as the
        stack(n_stacks).j = j;      % index for the current column
        stack(n_stacks).k_min = k;

        % Go through the column until we reach the end or find an exterior
        % point
        while (k <= N-1) && (grid(i,j,k).on == TRUE) 
           k = k + 1; 
        end
        stack(n_stacks).k_max = k - 1;

        % Add boundary terms for first and last points in current column
        stack(n_stacks).btf.x = x_min + hx * i;
        stack(n_stacks).btf.y = y_min + hy * j;
        stack(n_stacks).btf.z = zeta1(stack(n_stacks).btf.x, ...
            stack(n_stacks).btf.y);
        stack(n_stacks).btf.h_prime = z_min + hz * stack(n_stacks).k_min ...
            - stack(n_stacks).btf.z;    
        stack(n_stacks).btf.U = g3(stack(n_stacks).btf.x, stack(n_stacks).btf.y, ...
            stack(n_stacks).btf.z);
                % use initial conditions for first value of U

        % Add boundary terms for first and last points in current column
        stack(n_stacks).btl.x = x_min + hx * i;
        stack(n_stacks).btl.y = y_min + hy * j;
        stack(n_stacks).btl.z = zeta2(stack(n_stacks).btl.x, ...
            stack(n_stacks).btl.y);
        stack(n_stacks).btl.h_prime = stack(n_stacks).btl.z - ...
            (z_min + hz * stack(n_stacks).k_max);
        stack(n_stacks).btl.U = g3(stack(n_stacks).btl.x, stack(n_stacks).btl.y, ...
            stack(n_stacks).btl.z);
                % use initial conditions for first value of U
    end
end

% Define temporal grid.  Use a timestep of the same size as the larger of
% the two spactial grids. Go to time 1.
tau = max([hx hy hz]);
M = floor(1./tau) + 1;

% Allocate space for plots.
plot_data = zeros(N-1,N-1,9);

% Save initial state to plot for debugging
k = floor(N/2.);
for i = 1:N - 1
    for j = 1:N - 1
        if grid(i,j,k).on == 1
            plot_data(i,j,1) = grid(i,j,k).U; 
            plot_data(i,j,2) = u(grid(i,j,k).x, ...
                grid(i,j,k).y, 0);
            plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
        end
    end
end

% Create plots
X = linspace(x_min + hx, x_max - hx, N - 1);
Y = linspace(y_min + hy, y_max - hy, N - 1);

%Extra plots for debuging
figure
subplot(3,3,1)
surf(X,Y,plot_data(:,:,1))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution at initial condition')

subplot(3,3,2)
surf(X,Y,plot_data(:,:,2))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution at initial condition')

subplot(3,3,3)
surf(X,Y,plot_data(:,:,3));
colormap winter;
xlabel('x')
ylabel('y')
title('V at initial condition')

% subplot(3,3,4)
% surf(X,Y,plot_data(:,:,4))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Approximate solution after one timestep')
% 
% subplot(3,3,5)
% surf(X,Y,plot_data(:,:,5))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Exact solution after one timestep')
% 
% subplot(3,3,6)
% surf(X,Y,plot_data(:,:,6));
% colormap winter;
% xlabel('x')
% ylabel('y')
% title('Error after one timestep')
% 
% subplot(3,3,7)
% surf(X,Y,plot_data(:,:,7))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Approximate solution after final timestep')
% 
% subplot(3,3,8)
% surf(X,Y,plot_data(:,:,8))
% colormap winter
% xlabel('x')
% ylabel('y')
% title('Exact solution after final timestep')
% 
% subplot(3,3,9)
% surf(X,Y,plot_data(:,:,9));
% colormap winter;
% xlabel('x')
% ylabel('y')
% title('Error after final timestep')

end

