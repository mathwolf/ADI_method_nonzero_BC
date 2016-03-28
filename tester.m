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
               grid(i,j,k).r = 0;           % row
               grid(i,j,k).c = 0;           % col
               grid(i,j,k).s = 0;           % stack
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
            grid(i,j,k).r = n_rows;
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
        row(n_rows).btf.V = 0;
        
        row(n_rows).btl.y = y_min + hy * j;
        row(n_rows).btl.z = z_min + hz * k;
        
        row(n_rows).btl.x = eta2(row(n_rows).btl.y, row(n_rows).btl.z);
        row(n_rows).btl.h_prime = row(n_rows).btl.x ...
            -  ( x_min + hx * row(n_rows).i_max );   
        row(n_rows).btl.U = g3(row(n_rows).btl.x, row(n_rows).btl.y, ...
            row(n_rows).btl.z);
                % using initial conditions
        row(n_rows).btl.V = 0;
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
            grid(i,j,k).c = n_cols;
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
        col(n_cols).btf.V = 0;
        
        col(n_cols).btl.x = x_min + hx * i;
        col(n_cols).btl.z = z_min + hz * k;
        col(n_cols).btl.y = theta2(col(n_cols).btl.x, col(n_cols).btl.z);
        col(n_cols).btl.h_prime = col(n_cols).btl.y ...
            - ( y_min + hy * col(n_cols).j_max  );
        col(n_cols).btl.U = g3(col(n_cols).btl.x, col(n_cols).btl.y, ...
            col(n_cols).btl.z);
                % use initial conditions
        col(n_cols).btl.V = 0;
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
            grid(i,j,k).s = n_stacks;
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
        stack(n_stacks).btf.V = 0;
        
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
        stack(n_stacks).btl.V = 0;
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
            plot_data(i,j,2) = u3(grid(i,j,k).x, ...
                grid(i,j,k).y, grid(i,j,k).z, 0);
            plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
        end
    end
end

for m = 1:M
        % Update the RHS of step 1.  Here we need derivatives in all
    % three directions
    for i = 1:N-1
        for j = 1:N-1
           for k = 1:N-1
              if grid(i,j,k).on == TRUE
                 grid(i,j,k).V = grid(i,j,k).U;
                 % add the effect of the x derivative
                 % first consider the case where we are at the left
                 % boundary
                 r = grid(i,j,k).r;
                 if row(r).i_min == i
                     % we are next to the left boundary
                     if row(r).i_max == i
                         % we are also next to the right boundary
                         grid(i,j,k).V = grid(i,j,k).V + ...
                               (tau / (row(r).btf.h_prime + row(r).btl.h_prime) ) * ...
                               ( (row(r).btl.U - grid(i,j,k).U) / ...
                               row(r).btl.h_prime ...
                               - (grid(i,j,k).U - row(r).btf.U) / ...
                               row(r).btf.h_prime );
                     else
                         % left boundary but not right
                         grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (row(r).btf.h_prime + hx) ) * ...
                             ( (grid(i+1,j,k).U - grid(i,j,k).U) / hx ...
                             - (grid(i,j,k).U - row(r).btf.U) / ...
                             row(r).btf.h_prime );
                     end
                 else 
                     if row(r).i_max == i
                        % at the right boundary but not the left
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (hx + row(r).btl.h_prime) ) * ...
                             ( (row(r).btl.U - grid(i,j,k).U) / ...
                             row(r).btl.h_prime ...
                             - (grid(i,j,k).U - grid(i-1,j,k).U) / hx);
                     else
                         % point is not next to either boundary in x
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / 2*hx^2) * ...
                             (grid(i+1,j,k).U - 2*grid(i,j,k).U + ...
                             grid(i-1,j,k).U);
                     end
                 end
                 
                 % add the effect of the y derivative
                 % first consider the case where we are at the left
                 % boundary
                 c = grid(i,j,k).c;
                 if col(c).j_min == j
                     % we are next to the left boundary
                     if col(c).j_max == j
                         % we are also next to the right boundary
                         grid(i,j,k).V = grid(i,j,k).V + ...
                               (tau / (col(c).btf.h_prime + col(c).btl.h_prime) ) * ...
                               ( (col(c).btl.U - grid(i,j,k).U) / ...
                               col(c).btl.h_prime ...
                               - (grid(i,j,k).U - col(c).btf.U) / ...
                               col(c).btf.h_prime );
                     else
                         % left boundary but not right
                         grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (col(c).btf.h_prime + hy) ) * ...
                             ( (grid(i,j+1,k).U - grid(i,j,k).U) / hy ...
                             - (grid(i,j,k).U - col(c).btf.U) / ...
                             col(c).btf.h_prime );
                     end
                 else 
                     if col(c).j_max == j
                        % at the right boundary but not the left
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (hy + col(c).btl.h_prime) ) * ...
                             ( (col(c).btl.U - grid(i,j,k).U) / ...
                             col(c).btl.h_prime ...
                             - (grid(i,j,k).U - grid(i,j-1,k).U) / hy);
                     else
                         % point is not next to either boundary in x
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / 2*hy^2) * ...
                             (grid(i,j+1,k).U - 2*grid(i,j,k).U + ...
                             grid(i,j-1,k).U);
                     end
                 end
                 
                 % add the effect of the z derivative
                 % first consider the case where we are at the left
                 % boundary
                 s = grid(i,j,k).s;
                 if stack(s).k_min == k
                     % we are next to the left boundary
                     if stack(s).k_max == k
                         % we are also next to the right boundary
                         grid(i,j,k).V = grid(i,j,k).V + ...
                               (tau / (stack(s).btf.h_prime + stack(s).btl.h_prime) ) * ...
                               ( (stack(s).btl.U - grid(i,j,k).U) / ...
                               stack(s).btl.h_prime ...
                               - (grid(i,j,k).U - stack(s).btf.U) / ...
                               stack(s).btf.h_prime );
                     else
                         % left boundary but not right
                         grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (stack(s).btf.h_prime + hz) ) * ...
                             ( (grid(i,j,k+1).U - grid(i,j,k).U) / hz ...
                             - (grid(i,j,k).U - stack(s).btf.U) / ...
                             stack(s).btf.h_prime );
                     end
                 else 
                     if stack(s).k_max == k
                        % at the right boundary but not the left
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / (hz + stack(s).btl.h_prime) ) * ...
                             ( (stack(s).btl.U - grid(i,j,k).U) / ...
                             stack(s).btl.h_prime ...
                             - (grid(i,j,k).U - grid(i,j,k-1).U) / hz);
                     else
                         % point is not next to either boundary in x
                        grid(i,j,k).V = grid(i,j,k).V + ...
                             (tau / 2*hz^2) * ...
                             (grid(i,j,k+1).U - 2*grid(i,j,k).U + ...
                             grid(i,j,k-1).U);
                     end
                 end
                 
              end
           end
        end
    end
    
    % On the 1/3 timestep, we solve one matrix-vector equation for
    % each row
    for r = 1:n_rows
       j = row(r).j;
       k = row(r).k;
       i_min = row(r).i_min;
       i_max = row(r).i_max;
       length = i_max - i_min + 1;
       
       % Create vectors for rhs of equation
       b = zeros(length,1);
       for i = i_min:i_max
          b(i - i_min + 1) = grid(i,j,k).V; 
       end
              
       % Update the boundary values for the current row
       row(r).btf.V = g4(row(r).btf.x, row(r).btf.y, row(r).btf.z, tau*m);
       row(r).btl.V = g4(row(r).btl.x, row(r).btl.y, row(r).btl.z, tau*m);
       btf.W = row(r).btf.U + row(r).btf.V;
       btf.h_prime = row(r).btf.h_prime;
       btl.W = row(r).btl.U + row(r).btl.V;
       btl.h_prime = row(r).btl.h_prime;       
       
       % Add effect of boundary pts to rhs of equation
       % First, case where there is only one interior point in the row
       if length == 1
            b(1) = b(1) + ...
                (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * btf.W + ...
                (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * btl.W;
       % Case where A is is 2 x 2 or larger
       else
            b(1) = b(1) + (tau/((hx + btf.h_prime) * btf.h_prime) ) * btf.W;
            b(length) = b(length) + ...
                (tau/((hx + btl.h_prime) * btl.h_prime) ) * btl.W;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*hx^2)) * A;
       
       % Adjust first and last row of matrix since boundary points are 
       % unevenly spaced.  First, case when A is 1 x 1
       if length == 1
           A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);           
       % Case where A is 2 x 2 or larger
       else
            A(1,1) = tau / (hx * btf.h_prime);
            A(1,2) = - tau / (hx * (hx + btf.h_prime));
            A(length, length - 1) = - tau / (hx * (hx + btl.h_prime));
            A(length, length) = tau / (hx * btl.h_prime);
            A = speye(length) + A ;
       end
       
       % Solve the vector equation
       b = A\b;
       
       % Store the values of V for the current row
       for i = i_min:i_max
           grid(i,j,k).V = b(i - i_min + 1);
       end                   
    end
    
    % On the 2/3 timestep, we solve one matrix-vector equation for
    % each col
    for c = 1:n_cols
       i = col(c).i;
       k = col(c).k;
       j_min = col(c).j_min;
       j_max = col(c).j_max;
       length = j_max - j_min + 1;
       
       % Create vectors for rhs of equation
       b = zeros(length,1);
       for j = j_min:j_max
          b(j - j_min + 1) = grid(i,j,k).V; 
       end
                     
       % Update the boundary values for the current col
       col(c).btf.V = g4(col(c).btf.x, col(c).btf.y, col(c).btf.z, tau*m);
       col(c).btl.V = g4(col(c).btl.x, col(c).btl.y, col(c).btl.z, tau*m);
       btf.W = col(c).btf.V - col(c).btf.U;
       btf.h_prime = row(r).btf.h_prime;
       btl.W = row(r).btl.V - row(r).btl.U;
       btl.h_prime = row(r).btl.h_prime;       
       
       % Add effect of boundary pts to rhs of equation
       % First, case where there is only one interior point in the row
       if length == 1
            b(1) = b(1) + ...
                (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * btf.W + ...
                (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * btl.W;
       % Case where A is is 2 x 2 or larger
       else
            b(1) = b(1) + (tau/((hy + btf.h_prime) * btf.h_prime) ) * btf.W;
            b(length) = b(length) + ...
                (tau/((hy + btl.h_prime) * btl.h_prime) ) * btl.W;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*hy^2)) * A;
       
       % Adjust first and last row of matrix since boundary points are 
       % unevenly spaced.  First, case when A is 1 x 1
       if length == 1
           A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);           
       % Case where A is 2 x 2 or larger
       else
            A(1,1) = tau / (hy * btf.h_prime);
            A(1,2) = - tau / (hy * (hy + btf.h_prime));
            A(length, length - 1) = - tau / (hy * (hy + btl.h_prime));
            A(length, length) = tau / (hy * btl.h_prime);
            A = speye(length) + A ;
       end
       
       % Solve the vector equation
       b = A\b;
       
       % Store the values of V for the current col
       for j = j_min:j_max
           grid(i,j,k).V = b(j - j_min + 1);
       end                   
    end

    % Find the approximation to the PDE at the next whole timestep
    % Solve one matrix-vector equation for each stack   
    for s = 1:n_stacks
       i = stack(s).i;
       j = stack(s).j;
       k_min = stack(s).k_min;
       k_max = stack(s).k_max;
       length = k_max - k_min + 1;
       
       % Create vectors for rhs of equation
       b = zeros(length,1);
       for k = k_min:k_max
          b(k - k_min + 1) = grid(i,j,k).V; 
       end
                     
       % Update the boundary values for the current stack
       stack(s).btf.V = g4(stack(s).btf.x, stack(s).btf.y, stack(s).btf.z, tau*m);
       stack(s).btl.V = g4(stack(s).btl.x, stack(s).btl.y, stack(s).btl.z, tau*m);
       btf.W = col(c).btf.V - col(c).btf.U;
       btf.h_prime = row(r).btf.h_prime;
       btl.W = row(r).btl.V - row(r).btl.U;
       btl.h_prime = row(r).btl.h_prime;       
       
       % Add effect of boundary pts to rhs of equation
       % First, case where there is only one interior point in the row
       if length == 1
            b(1) = b(1) + ...
                (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * btf.W + ...
                (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * btl.W;
       % Case where A is is 2 x 2 or larger
       else
            b(1) = b(1) + (tau/((hz + btf.h_prime) * btf.h_prime) ) * btf.W;
            b(length) = b(length) + ...
                (tau/((hz + btl.h_prime) * btl.h_prime) ) * btl.W;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*hz^2)) * A;
       
       % Adjust first and last row of matrix since boundary points are 
       % unevenly spaced.  First, case when A is 1 x 1
       if length == 1
           A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);           
       % Case where A is 2 x 2 or larger
       else
            A(1,1) = tau / (hz * btf.h_prime);
            A(1,2) = - tau / (hz * (hz + btf.h_prime));
            A(length, length - 1) = - tau / (hz * (hz + btl.h_prime));
            A(length, length) = tau / (hz * btl.h_prime);
            A = speye(length) + A ;
       end
       
       % Solve the vector equation
       b = A\b;
       
       % Store the values of U for the current stack
       for k = k_min:k_max
           grid(i,j,k).U = b(k - k_min + 1);
       end                   
    end
    
    

end

% Create plots
X = linspace(x_min + hx, x_max - hx, N - 1);
Y = linspace(y_min + hy, y_max - hy, N - 1);

%Extra plots for debuging
figure
subplot(1,3,1)
surf(X,Y,plot_data(:,:,1))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution at initial condition')

subplot(1,3,2)
surf(X,Y,plot_data(:,:,2))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution at initial condition')

subplot(1,3,3)
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
