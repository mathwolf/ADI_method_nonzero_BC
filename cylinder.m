function cylinder()

% Matlab R2013a

% ADI method for the solution of the parabolic PDE with Dirichlet boundary
% conditions.  In this method the domain of the problem is an arbitrary,
% possibly nonrectangular region.

TRUE = 1;
FALSE = 0;

% Constants used to switch between different test functions.
EXPONENT_0 = 0;
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLYNOMIAL = 4;
EXPONENT_3 = 5;
global test_solution
test_solution = EXPONENT_1;

% Constants used to switch between different test domains.
SPHERE = 1;
ELLIPSOID = 2;
OCTAHEDRON = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
CUBE = 7;

CIRCLE_CYLINDER = 8;
ELLIPSE_CYLINDER = 9;

global domain
domain = ELLIPSE_CYLINDER;

% Description of spatial grid.  We divide both the x and y dimensions of
% the problem into the same number of gridpoints.  First, identify the min 
% and max values in both directions.  The actual domain of the PDE will be 
% a subset of the big rectangle defined by these values.

% For now, we are just defaulting to a cube of side length 2.
if domain == ELLIPSOID
   x_min = -1;
   x_max = 1;
   y_min = -0.5;
   y_max = 0.5;
   z_min = -0.25;
   z_max = 0.25;
elseif domain == CUBE
   x_min = 0;
   x_max = 1;
   y_min = 0;
   y_max = 1;
   z_min = 0;
   z_max = 1;
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
elseif domain == CIRCLE_CYLINDER
    x_min = -1;
    x_max = 1;
    y_min = -1;
    y_max = 1;
    z_min = 0;
    z_max = 1;
elseif domain == ELLIPSE_CYLINDER
    x_min = -1;
    x_max = 1;
    y_min = -0.5;
    y_max = 0.5;
    z_min = 0;
    z_max = 1;
else
    % all the other test domains have the same spatial range
   x_min = -1;
   x_max = 1;
   y_min = -1;
   y_max = 1;
   z_min = -1;
   z_max = 1;
end

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the size by 
% half.

% Check three different test functions
for q = 2:4
    if q == 4
        q = 5;
    end
    test_solution = q;
    disp('Test function');
    disp(q);

for p = 1:5
    
    % Use a spatial grid of 0.4 times 2 to the power p-1
    N = 5 * 2^(p-1);
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
                % Include possible round-off error when checking interior points
                if (phi1(x) + 10^(-10) < y) && ...
                   (y < phi2(x) - 10^(-10)) && ...
                   (zeta1(x,y) + 10^(-10) < z) && ...
                   (z < zeta2(x,y) - 10^(-10)) && ...
                   (eta1(y,z) + 10^(-10) < x) && ...
                   (x < eta2(y,z) - 10^(-10)) && ...
                   (theta1(x,z) + 10^(-10) < y) && ...
                   (y < theta2(x,z) - 10^(-10))
                   
                   % The point is an interio point
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

    % Create a data structure that describes the rows of the grid.
    n_rows = 0;
    for k = 1:N-1
        for j = 1:N-1
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

            % Add boundary terms for first and last points in current stack
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
%      tau = max([hx hy hz]);
%      M = floor(1./tau) + 1;
    tau = 1. / N;
    M = N;
    
%    For debugging
%     for k = 1:N-1
%        for j = 1:N-1
%           for i = 1:N-1
%              disp(i);
%              disp(j);
%              disp(k);
%              disp(grid(i,j,k)); 
%           end
%        end
%     end
    
    % Step through the scheme
    for m = 1:M
        % Update the RHS of step 1.  Here we need derivatives in all
        % three directions
        for i = 1:N-1
            for j = 1:N-1
               for k = 1:N-1
                  if grid(i,j,k).on == TRUE
                     % add the effect of the x derivative
                     % first consider the case where we are at the left
                     % boundary
                     r = grid(i,j,k).r;
                     if row(r).i_min == i
                         % we are next to the left boundary
                         if row(r).i_max == i
                             % we are also next to the right boundary
                             grid(i,j,k).V = grid(i,j,k).U + ...
                                   (tau / (row(r).btf.h_prime + row(r).btl.h_prime) ) * ...
                                   ( (row(r).btl.U - grid(i,j,k).U) / ...
                                   row(r).btl.h_prime ...
                                   - (grid(i,j,k).U - row(r).btf.U) / ...
                                   row(r).btf.h_prime );
                         else
                             % left boundary but not right
                             grid(i,j,k).V = grid(i,j,k).U + ...
                                 (tau / (row(r).btf.h_prime + hx) ) * ...
                                 ( (grid(i+1,j,k).U - grid(i,j,k).U) / hx ...
                                 - (grid(i,j,k).U - row(r).btf.U) / ...
                                 row(r).btf.h_prime );
                         end
                     else 
                         if row(r).i_max == i
                            % at the right boundary but not the left
                            grid(i,j,k).V = grid(i,j,k).U + ...
                                 (tau / (hx + row(r).btl.h_prime) ) * ...
                                 ( (row(r).btl.U - grid(i,j,k).U) / ...
                                 row(r).btl.h_prime ...
                                 - (grid(i,j,k).U - grid(i-1,j,k).U) / hx);
                         else
                             % point is not next to either boundary in x
                            grid(i,j,k).V = grid(i,j,k).U + ...
                                 (tau / (2*hx^2)) * ...
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
                                   (2 * tau / (col(c).btf.h_prime + col(c).btl.h_prime) ) * ...
                                   ( (col(c).btl.U - grid(i,j,k).U) / ...
                                   col(c).btl.h_prime ...
                                   - (grid(i,j,k).U - col(c).btf.U) / ...
                                   col(c).btf.h_prime );
                         else
                             % left boundary but not right
                             grid(i,j,k).V = grid(i,j,k).V + ...
                                 (2 * tau / (col(c).btf.h_prime + hy) ) * ...
                                 ( (grid(i,j+1,k).U - grid(i,j,k).U) / hy ...
                                 - (grid(i,j,k).U - col(c).btf.U) / ...
                                 col(c).btf.h_prime );
                         end
                     else 
                         if col(c).j_max == j
                            % at the right boundary but not the left
                            grid(i,j,k).V = grid(i,j,k).V + ...
                                 (2 * tau / (hy + col(c).btl.h_prime) ) * ...
                                 ( (col(c).btl.U - grid(i,j,k).U) / ...
                                 col(c).btl.h_prime ...
                                 - (grid(i,j,k).U - grid(i,j-1,k).U) / hy);
                         else
                             % point is not next to either boundary in x
                            grid(i,j,k).V = grid(i,j,k).V + ...
                                 (tau / hy^2) * ...
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
                                   (2 * tau / (stack(s).btf.h_prime + stack(s).btl.h_prime) ) * ...
                                   ( (stack(s).btl.U - grid(i,j,k).U) / ...
                                   stack(s).btl.h_prime ...
                                   - (grid(i,j,k).U - stack(s).btf.U) / ...
                                   stack(s).btf.h_prime );
                         else
                             % left boundary but not right
                             grid(i,j,k).V = grid(i,j,k).V + ...
                                 (2 * tau / (stack(s).btf.h_prime + hz) ) * ...
                                 ( (grid(i,j,k+1).U - grid(i,j,k).U) / hz ...
                                 - (grid(i,j,k).U - stack(s).btf.U) / ...
                                 stack(s).btf.h_prime );
                         end
                     else 
                         if stack(s).k_max == k
                            % at the right boundary but not the left
                            grid(i,j,k).V = grid(i,j,k).V + ...
                                 (2 * tau / (hz + stack(s).btl.h_prime) ) * ...
                                 ( (stack(s).btl.U - grid(i,j,k).U) / ...
                                 stack(s).btl.h_prime ...
                                 - (grid(i,j,k).U - grid(i,j,k-1).U) / hz);
                         else
                             % point is not next to either boundary in x
                            grid(i,j,k).V = grid(i,j,k).V + ...
                                 (tau / hz^2) * ...
                                 (grid(i,j,k+1).U - 2*grid(i,j,k).U + ...
                                 grid(i,j,k-1).U);
                         end
                     end
                  end
                    end
                end
            end

            for i = 1:N-1
                for j = 1:N-1
                    for k = 1:N-1
                        if (grid(i,j,k).on == TRUE)
                     % Add the effect of the forcing term
                     x = grid(i,j,k).x;
                     y = grid(i,j,k).y;
                     z = grid(i,j,k).z;
                     grid(i,j,k).V = grid(i,j,k).V + ...
                         (tau/2)*(f3(x,y,z,tau*(m-1)) + ...
                         f3(x,y,z,tau*m));
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

           % Update the boundary values for the current row,
           % using partial perturbation in the z-direction
           row(r).btf.V = g4(row(r).btf.x, row(r).btf.y, row(r).btf.z, tau*m) ...
                - (tau/(2*hz^2)) * (g4(row(r).btf.x, row(r).btf.y, row(r).btf.z - hz, tau*m) ...
                    - 2*g4(row(r).btf.x, row(r).btf.y, row(r).btf.z, tau*m) ...
                    + g4(row(r).btf.x, row(r).btf.y, row(r).btf.z + hz, tau*m)) ...
                + (tau/(2*hz^2)) * (g4(row(r).btf.x, row(r).btf.y, row(r).btf.z - hz, tau*(m-1)) ...
                    - 2*g4(row(r).btf.x, row(r).btf.y, row(r).btf.z, tau*(m-1)) ...
                    + g4(row(r).btf.x, row(r).btf.y, row(r).btf.z + hz, tau*(m-1)));
           row(r).btl.V = g4(row(r).btl.x, row(r).btl.y, row(r).btl.z, tau*m) ...
                - (tau/(2*hz^2)) * (g4(row(r).btl.x, row(r).btl.y, row(r).btl.z - hz, tau*m) ...
                    - 2*g4(row(r).btl.x, row(r).btl.y, row(r).btl.z, tau*m) ...
                    + g4(row(r).btl.x, row(r).btl.y, row(r).btl.z + hz, tau*m)) ...
                + (tau/(2*hz^2)) * (g4(row(r).btl.x, row(r).btl.y, row(r).btl.z - hz, tau*(m-1)) ...
                    - 2*g4(row(r).btl.x, row(r).btl.y, row(r).btl.z, tau*(m-1)) ...
                    + g4(row(r).btl.x, row(r).btl.y, row(r).btl.z + hz, tau*(m-1)));
           btf.h_prime = row(r).btf.h_prime;
           btl.h_prime = row(r).btl.h_prime;       

           % Add effect of boundary pts to rhs of equation
           % First, case where there is only one interior point in the row
           if length == 1
                b(1) = b(1) + ...
                    (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * row(r).btf.V + ...
                    (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * row(r).btl.V;
           % Case where A is is 2 x 2 or larger
           else
                b(1) = b(1) + (tau/((hx + btf.h_prime) * btf.h_prime) ) * row(r).btf.V;
                b(length) = b(length) + ...
                    (tau/((hx + btl.h_prime) * btl.h_prime) ) * row(r).btl.V;
           end

           % Set up tridiagonal matrix on lhs of equation
           A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                            [-1,0,1], length, length);
           A = (tau/(2*hx^2)) * A;

           % Adjust first and last row of matrix since boundary points are 
           % unevenly spaced.  First, case when A is 1 x 1
           if length == 1
               A = 1 + tau / (btf.h_prime * btl.h_prime);           
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
           % Store new boundary values for current row
           row(r).btf.U = row(r).btf.V;
           row(r).btl.U = row(r).btl.V;
        end

        % On the 2/3 timestep, we solve one matrix-vector equation for
        % each col
        for c = 1:n_cols
           i = col(c).i;
           k = col(c).k;
           j_min = col(c).j_min;
           j_max = col(c).j_max;
           length = j_max - j_min + 1;

           % Set up tridiagonal matrix for lhs of equation
           tridiagonal = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                            [-1,0,1], length, length);
           tridiagonal = (tau/(2*hy^2)) * tridiagonal;

           % Create vectors for rhs of equation
           b = zeros(length,1);
           d = zeros(length,1);
           for j = j_min:j_max
              b(j - j_min + 1) = grid(i,j,k).V; 
              d(j - j_min + 1) = grid(i,j,k).U;
           end

           % Add effects of y-derivative
           b = b + tridiagonal*d;

           % Update the boundary values for the current col
           col(c).btf.V = g4(col(c).btf.x, col(c).btf.y, col(c).btf.z, tau*m) ...
                - (tau/(2*hz^2)) * (g4(col(c).btf.x, col(c).btf.y, col(c).btf.z - hz, tau*m) ...
                    - 2*g4(col(c).btf.x, col(c).btf.y, col(c).btf.z, tau*m) ...
                    + g4(col(c).btf.x, col(c).btf.y, col(c).btf.z + hz, tau*m)) ...
                + (tau/(2*hz^2)) * (g4(col(c).btf.x, col(c).btf.y, col(c).btf.z - hz, tau*(m-1)) ...
                    - 2*g4(col(c).btf.x, col(c).btf.y, col(c).btf.z, tau*(m-1)) ...
                    + g4(col(c).btf.x, col(c).btf.y, col(c).btf.z + hz, tau*(m-1)));
           col(c).btl.V = g4(col(c).btl.x, col(c).btl.y, col(c).btl.z, tau*m) ...
                - (tau/(2*hz^2)) * (g4(col(c).btl.x, col(c).btl.y, col(c).btl.z - hz, tau*m) ...
                    - 2*g4(col(c).btl.x, col(c).btl.y, col(c).btl.z, tau*m) ...
                    + g4(col(c).btl.x, col(c).btl.y, col(c).btl.z + hz, tau*m)) ...
                + (tau/(2*hz^2)) * (g4(col(c).btl.x, col(c).btl.y, col(c).btl.z - hz, tau*(m-1)) ...
                    - 2*g4(col(c).btl.x, col(c).btl.y, col(c).btl.z, tau*(m-1)) ...
                    + g4(col(c).btl.x, col(c).btl.y, col(c).btl.z + hz, tau*(m-1)));
           btf.h_prime = col(c).btf.h_prime;
           btl.h_prime = col(c).btl.h_prime;       

           % Add effect of boundary pts to rhs of equation
           % First, case where there is only one interior point in the row
           if length == 1
                b(1) = grid(i,j_min,k).V - tau/(btf.h_prime + btl.h_prime) *...
                    ( (col(c).btf.U - grid(i,j_min,k).U)/btf.h_prime - ...
                    (grid(i,j_min,k).U - col(c).btl.U)/btl.h_prime) + ...
                    tau/(btf.h_prime + btl.h_prime) * ...
                    ( col(c).btf.V/btf.h_prime + col(c).btl.V/btl.h_prime );
           % Case where A is is 2 x 2 or larger
           else
                b(1) = grid(i,j_min,k).V - tau/(hy + btf.h_prime) * ...
                    ( (col(c).btf.U - grid(i,j_min,k).U) / btf.h_prime - ...
                    (grid(i,j_min,k).U - grid(i,j_min+1,k).U) / hy ) + ...
                    tau/((hy + btf.h_prime)*btf.h_prime) * col(c).btf.V;
                b(length) = grid(i,j_max,k).V - tau/(hy + btl.h_prime) * ...
                    ( (col(c).btl.U - grid(i,j_max,k).U) / btl.h_prime - ...
                    (grid(i,j_max,k).U - grid(i,j_max-1,k).U) / hy ) + ...
                    tau/((hy + btl.h_prime)*btl.h_prime) * col(c).btl.V;
           end

           % Adjust first and last row of matrix since boundary points are 
           % unevenly spaced.  First, case when A is 1 x 1
           if length == 1
               A = 1 + tau / (btf.h_prime * btl.h_prime);           
           % Case where A is 2 x 2 or larger
           else
                A = tridiagonal;
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
           % Store boundary values for the current col
           col(c).btf.U = col(c).btf.V;
           col(c).btl.U = col(c).btl.V;
        end

        % Find the approximation to the PDE at the next whole timestep
        % Solve one matrix-vector equation for each stack   
        for s = 1:n_stacks
           i = stack(s).i;
           j = stack(s).j;
           k_min = stack(s).k_min;
           k_max = stack(s).k_max;
           length = k_max - k_min + 1;

           % Set up tridiagonal matrix on lhs of equation
           tridiagonal = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                            [-1,0,1], length, length);
           tridiagonal = (tau/(2*hz^2)) * tridiagonal;

           % Create vectors for rhs of equation
           b = zeros(length,1);
           d = zeros(length,1);
           for k = k_min:k_max
              b(k - k_min + 1) = grid(i,j,k).V; 
              d(k - k_min + 1) = grid(i,j,k).U;
           end

           % Add effects of the z-derivative
           b = b + tridiagonal * d;

           % Update the boundary values for the current stack
           stack(s).btf.V = g4(stack(s).btf.x, stack(s).btf.y, stack(s).btf.z, tau*m);
           stack(s).btl.V = g4(stack(s).btl.x, stack(s).btl.y, stack(s).btl.z, tau*m);
           btf.h_prime = stack(s).btf.h_prime;
           btl.h_prime = stack(s).btl.h_prime;       

           % Add effect of boundary pts to rhs of equation
           % First, case where there is only one interior point in the row
           if length == 1
                b(1) = grid(i,j,k_min).V - tau/(btf.h_prime + btl.h_prime) * ...
                    ( (stack(s).btf.U - grid(i,j,k_min).U)/btf.h_prime - ...
                    (grid(i,j,k_min).U - stack(s).btl.U)/btl.h_prime) + ...
                    tau/(btf.h_prime + btl.h_prime) * ...
                    ( stack(s).btf.V/btf.h_prime + stack(s).btl.V/btl.h_prime );
           % Case where A is is 2 x 2 or larger
           else
                b(1) = grid(i,j,k_min).V - tau/(hz + btf.h_prime) * ...
                    ( (stack(s).btf.U - grid(i,j,k_min).U) / btf.h_prime - ...
                    (grid(i,j,k_min).U - grid(i,j,k_min+1).U) / hz ) + ...
                    tau/((hz + btf.h_prime)*btf.h_prime) * stack(s).btf.V;
                b(length) = grid(i,j,k_max).V - tau/(hz + btl.h_prime) * ...
                    ( (stack(s).btl.U - grid(i,j,k_max).U) / btl.h_prime - ...
                    (grid(i,j,k_max).U - grid(i,j,k_max-1).U) / hz ) + ...
                    tau/((hz + btl.h_prime)*btl.h_prime) * stack(s).btl.V;
           end

           % Adjust first and last row of matrix since boundary points are 
           % unevenly spaced.  First, case when A is 1 x 1
           if length == 1
               A = 1 + tau / (btf.h_prime * btl.h_prime);           
           % Case where A is 2 x 2 or larger
           else
                A = tridiagonal;
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
           % Store boundary values for the current stack
           stack(s).btf.U = stack(s).btf.V;
           stack(s).btl.U = stack(s).btl.V;
        end
    end

    % Evaluate error at interior gridpoints at final timestep
    max_nodal_error = 0;
    nodal_error_squared = 0;    
    root_mean_sq_error = 0;
    debug_i = 0;
    debug_j = 0;
    debug_k = 0;

    for i = 1:N-1
       for j = 1:N-1
           for k = 1:N-1
               if grid(i,j,k).on == TRUE
                   current_error = abs(grid(i,j,k).U - ...
                       u3(grid(i,j,k).x, grid(i,j,k).y, grid(i,j,k).z, tau*m));
                   if current_error > max_nodal_error
                      max_nodal_error = current_error; 
                   end
                   nodal_error_squared = nodal_error_squared + current_error^2;
                   debug_i = i;
                   debug_j = j;
                   debug_k = k;
               end
           end
       end
    end
    root_mean_sq_error = sqrt(hx*hy*hz * nodal_error_squared);

%     disp(debug_i);
%     disp(debug_j);
%     disp(debug_k);
%     disp(99999);
    % Write data to table
    table_data(p,1) = tau;
    table_data(p,2) = max_nodal_error;
    if p ~= 1
       % estimate order of convergence for max err
       table_data(p,3) = log(table_data(p-1,2) / table_data(p,2)) / ...
           log(table_data(p-1,1) / table_data(p,1)); 
    end   
    table_data(p,4) = tau;
    table_data(p,5) = root_mean_sq_error;
    if p ~= 1
       % estimate order of convergence for rms err
       table_data(p,6) = log(table_data(p-1,5) / table_data(p,5)) / ...
           log(table_data(p-1,4) / table_data(p,4));
    end
    
end

% Display table of error information
disp('    h                  max error          order of conv');
disp(table_data(:,1:3));
disp('    h                  rms error          order of conv');
        disp(table_data(:,4:6));
end
end