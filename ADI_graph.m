function ADI_graph()

% Matlab R2013a

% ADI method for the solution of the parabolic PDE with Dirichlet boundary
% conditions.  In this method the domain of the problem is an arbitrary,
% possibly nonrectangular region.

TRUE = 1;
FALSE = 0;

% Constants used to switch between different test functions.
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLYNOMIAL = 4;
global test_solution
test_solution = TRIG;

% Constants used to switch between different test domains.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;
global domain
domain = CIRCLE;

% Description of spatial grid.  At this point, the method is set up to use 
% a grid with the same uniform spacing in x and y.  The grid goes
% from -1 to +1.  The actual domain of the PDE will be a subset of these
% gridpoints.
N= 41;
h = 1./20.;
x_min = -1.;
y_min = -1.;

% Define the gridpoints that are part of the interior of the domain of the
% PDE.  Check each gridpoint in the array to see if it is an interior 
% point.  Examine one column at a time beginning with the leftmost column
% of the grid.
for i = 1:N
   % Check each point in the column starting with the first row
   for j = 1:N
       x = x_min + h * (i-1);
       y = y_min + h * (j-1);
       if (phi1(x) < y) && (y < phi2(x)) && ...
                (psi1(y) < x) && (x < psi2(y))
           % The point is an interior point
           grid(i,j).on = TRUE; 
           grid(i,j).x = x;
           grid(i,j).y = y;
           grid(i,j).U = g1(x,y);   % use initial condition to fill in the
                                    % value of the approximation
           grid(i,j).V = 0.;
       else
           % The point is an exterior point
           grid(i,j).on = FALSE;
       end
   end
end

% Create a data structure that describes the rows of the grid.
n_rows = 0;
for j = 1:N
    % Go along the row until we reach the end or find an interior point.
    i = 1;
    while (i <= N) && (grid(i,j).on == FALSE)
       i = i + 1; 
    end

    if i > N
       % We are at the end of an empty row, go to the next row
       continue;
    end
    
    % We are at the first interior point in a new row
    n_rows = n_rows + 1;    
    row(n_rows).j = j;      % assign the next number in sequence as the
                            % index for the current row
    row(n_rows).i_min = i;
    
    % Go along the row until we reach the end or or find an exterior point
    while (i <= N) && (grid(i,j).on == TRUE)
       i = i + 1; 
    end
    row(n_rows).i_max = i - 1;
    
    % Add boundary terms for first and last point in current row
    row(n_rows).btf.y = y_min + h * (j-1);
    row(n_rows).btf.x = psi1(row(n_rows).btf.y);
    row(n_rows).btf.h_prime = x_min + h * (row(n_rows).i_min - 1) ...
        - row(n_rows).btf.x;
    row(n_rows).btf.U = g1(row(n_rows).btf.x, row(n_rows).btf.y);
            % using initial conditions for first value of U
    
    row(n_rows).btl.y = y_min + h * (j-1);
    row(n_rows).btl.x = psi2(row(n_rows).btl.y);
    row(n_rows).btl.h_prime = row(n_rows).btl.x ...
        - ( x_min + h * (row(n_rows).i_max - 1) );   
    row(n_rows).btl.U = g1(row(n_rows).btl.x, row(n_rows).btl.y);
            % using initial conditions
end

% Create a data structure that describes the columns of the grid.
n_cols = 0;
for i = 1:N
    % Go through col until we reach the end or find an interior point
    j = 1;
    while (j <= N) && (grid(i,j).on == FALSE) 
       j = j + 1; 
    end

    if j > N
       % We are at the end of an empty col, go to the next col
       continue;
    end
    
    % We are at the first interior point in a new col
    n_cols = n_cols + 1;
    col(n_cols).i = i;      % assign the next number in sequence as the
                            % index for the current column
    col(n_cols).j_min = j;
    
    % Go through the column until we reach the end or find an exterior
    % point
    while (j <= N) && (grid(i,j).on == TRUE) 
       j = j + 1; 
    end
    col(n_cols).j_max = j - 1;
    
    % Add boundary terms for first and last points in current column
    col(n_cols).btf.x = x_min + h * (i-1);
    col(n_cols).btf.y = phi1(col(n_cols).btf.x);
    col(n_cols).btf.h_prime = y_min + h * (col(n_cols).j_min - 1) ...
        - col(n_cols).btf.y;    
    col(n_cols).btf.U = g1(col(n_cols).btf.x, col(n_cols).btf.y);
            % use initial conditions for first value of U
            
    col(n_cols).btl.x = x_min + h * (i-1);
    col(n_cols).btl.y = phi2(col(n_cols).btl.x);
    col(n_cols).btl.h_prime = col(n_cols).btl.y ...
        - ( y_min + h * (col(n_cols).j_max - 1) );
    col(n_cols).btl.U = g1(col(n_cols).btl.x, col(n_cols).btl.y);
            % use initial conditions
end

% Define temporal grid.  Use a timestep of the same size as the space grid.
% Go to time 1.
tau = h;
M = 1./tau;

% Allocate space for plots.
plot_data = zeros(N,N,9);

% Save initial state to plot for debugging
for i = 1:N
    for j = 1:N
        if grid(i,j).on == 1
            plot_data(i,j,1) = grid(i,j).U; 
            plot_data(i,j,2) = u(grid(i,j).x, ...
                grid(i,j).y, 0);
            plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
        end
    end
end


% Evaluate each timestep.
for m = 1:M
    
    % Calculate the value of the function V at every interior gridpoint.
    for c = 1:n_cols
       i = col(c).i;
       j = col(c).j_min;
       
       % First, consider special case where the column holds only 1
       % interior point.
       if j == col(c).j_max
           grid(i,j).V = grid(i,j).U + ...
               (tau / (col(c).btf.h_prime + col(c).btl.h_prime) ) * ...
               ( (col(c).btl.U - grid(i,j).U) / ...
               col(c).btl.h_prime ...
               - (grid(i,j).U - col(c).btf.U) / ...
               col(c).btf.h_prime );
           
       % Otherwise the col holds two or more interior points
       else
            grid(i,j).V = grid(i,j).U + ...
                (tau / (col(c).btf.h_prime + h) ) * ...
                ( (grid(i,j+1).U - grid(i,j).U) / h ...
                - (grid(i,j).U - col(c).btf.U) / ...
                col(c).btf.h_prime );
            j = j + 1;
            while (j < col(c).j_max)
                grid(i,j).V = grid(i,j).U + ...
                    (tau / (2*h^2)) * (grid(i,j+1).U ...
                    - 2 * grid(i,j).U ...
                    + grid(i,j-1).U);
                j = j + 1;
            end
            grid(i,j).V = grid(i,j).U + ...
                (tau / (h + col(c).btl.h_prime) ) * ...
                ( (col(c).btl.U - grid(i,j).U) / ...
                col(c).btl.h_prime ...
                - (grid(i,j).U - grid(i,j-1).U) / h );                
       end
    end
    
    % Plot V for debuggging
    if m == 1
       for i = 1:N
          for j = 1:N
              if grid(i,j).on == TRUE
                 plot_data(i,j,3) = grid(i,j).V; 
              end
          end
       end
    end
    
    
    % Find the approximation U at the nest half-timestep
    % Here we solve a matrix-vector equation once for each row
    for r = 1:n_rows
       j = row(r).j;
       i_min = row(r).i_min;
       i_max = row(r).i_max;
       length = i_max - i_min + 1;
       
       % Create vectors for rhs of equation
       V = zeros(length,1);
       load_vector = zeros(length,1);
       for i = i_min:i_max
          V(i - i_min + 1) = grid(i,j).V; 
          load_vector(i - i_min + 1) = (0.5) * ...
              ( f(grid(i,j).x, grid(i,j).y, tau*(m-1)) + ...
              f(grid(i,j).x, grid(i,j).y, tau*m) );
       end
              
       % rhs of vector equation
       b = V + (tau/2) * load_vector;  
       
       % Update the boundary values for the current row
       row(r).btf.U = g2(row(r).btf.x, row(r).btf.y, tau*(m-0.5));
       row(r).btl.U = g2(row(r).btl.x, row(r).btl.y, tau*(m-0.5));
       btf.U = row(r).btf.U;
       btf.h_prime = row(r).btf.h_prime;
       btl.U = row(r).btl.U;
       btl.h_prime = row(r).btl.h_prime;       
       
       % Add effect of boundary pts to rhs of equation
       % First, case where A is 1 x 1
       if length == 1
            b(1) = b(1) + ...
                (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * btf.U + ...
                (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * btl.U;
       % Case where A is is 2 x 2 or larger
       else
            b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*h^2)) * A;
       
       % Adjust first and last row of matrix since boundary points are 
       % unevenly spaced.  First, case when A is 1 x 1
       if length == 1
           A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);           
       % Case where A is 2 x 2 or larger
       else
            A(1,1) = tau / (h * btf.h_prime);
            A(1,2) = - tau / (h * (h + btf.h_prime));
            A(length, length - 1) = - tau / (h * (h + btl.h_prime));
            A(length, length) = tau / (h * btl.h_prime);
            A = speye(length) + A ;
       end
       
       % Solve the vector equation
       x = A\b;
       
       % Store the values of U for the current row
       for i = i_min:i_max
           grid(i,j).U = x(i - i_min + 1);
       end                   
    end
    
    
    % Plot half-timestep solution for debugging
    if m == 1
       for i = 1:N
          for j = 1:N
              if grid(i,j).on == 1
                 plot_data(i,j,4) = grid(i,j).U; 
                 plot_data(i,j,5) = u(grid(i,j).x, ...
                     grid(i,j).y, tau*(m-0.5));
                 plot_data(i,j,6) = ...
                     plot_data(i,j,4) - plot_data(i,j,5);
              end
          end
       end
    end
    
    
    % Find the approximation U for the next whole timestep.
    % Here we solve a matrix/vector equation once for each column
    for c = 1:n_cols
        i = col(c).i;
        j_min = col(c).j_min;
        j_max = col(c).j_max;
        length = col(c).j_max - col(c).j_min + 1;
        
        % Create vectors for rhs of equation
        U = zeros(length,1);
        V = zeros(length,1);
        for j = j_min:j_max
            U(j - j_min + 1) = grid(i,j).U;
            V(j - j_min + 1) = grid(i,j).V;
        end
        
        % rhs of vector equation
        b = 2*U - V;
        
       % Update the boundary values for the current column
       col(c).btf.U = g2(col(c).btf.x, col(c).btf.y, tau*m);
       col(c).btl.U = g2(col(c).btl.x, col(c).btl.y, tau*m);
       btf.U = col(c).btf.U;
       btf.h_prime = col(c).btf.h_prime;
       btl.U = col(c).btl.U;
       btl.h_prime = col(c).btl.h_prime;
       
       % Add effect of boundary pts to rhs of equation
       % First, case where A is 1 x 1
       if length == 1
            b(1) = b(1) + ...
                (tau/((btf.h_prime + btl.h_prime) * btf.h_prime) ) * btf.U + ...
                (tau/((btf.h_prime + btl.h_prime) * btl.h_prime) ) * btl.U;
       % Case where A is is 2 x 2 or larger
       else
            b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*h^2)) * A;
       % Adjust first and last row since boundary points are unevenly
       % spaced.  First consider case when A is 1 x 1
       if length == 1
            A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);
       % Then case when A is 2 x 2 or larger
       else
            A(1,1) = tau / (h * btf.h_prime);
            A(1,2) = - tau / (h * (h + btf.h_prime));
            A(length, length - 1) = - tau / (h * (h + btl.h_prime));
            A(length, length) = tau / (h * btl.h_prime);       
            A = speye(length) + A ;
       end
       
       % Solve the vector equation
       x = A\b;

       % Store the values of U for the current column
       for j = j_min:j_max
           grid(i,j).U = x(j - j_min + 1);
       end
    end
    
    % Store data for plotting on selected steps
    if m == 1
       for i = 1:N
          for j = 1:N
              if grid(i,j).on == 1
                 plot_data(i,j,4) = grid(i,j).U; 
                 plot_data(i,j,5) = u(grid(i,j).x, ...
                     grid(i,j).y, tau*m);
                 plot_data(i,j,6) = ...
                     plot_data(i,j,4) - plot_data(i,j,5);
              end
          end
       end
       
    elseif m == M
        for i = 1:N
          for j = 1:N
              if grid(i,j).on == 1
                 plot_data(i,j,7) = grid(i,j).U; 
                 plot_data(i,j,8) = u(grid(i,j).x, ...
                     grid(i,j).y, tau*m);
                 plot_data(i,j,9) = ...
                     abs(plot_data(i,j,7) - plot_data(i,j,8));
              end
          end
       
        end

    end
end
% Create plots
X = linspace(-1., 1., N);
Y = linspace(-1., 1., N);


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

subplot(3,3,4)
surf(X,Y,plot_data(:,:,4))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after one half timestep')

subplot(3,3,5)
surf(X,Y,plot_data(:,:,5))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after one half timestep')

subplot(3,3,6)
surf(X,Y,plot_data(:,:,6));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after one half timestep')

subplot(3,3,7)
surf(X,Y,plot_data(:,:,7))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after one timestep')

subplot(3,3,8)
surf(X,Y,plot_data(:,:,8))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after one timestep')

subplot(3,3,9)
surf(X,Y,plot_data(:,:,9));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after one timestep')

end