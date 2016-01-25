function ADI_graph()

% Matlab R2013a

% ADI method for the solution of the parabolic PDE with Dirichlet boundary
% conditions.  
% Moving from 2D to 3D.  Start with 0 Dirichlet boundary conditions,
% a square domain, and the simple exponential test function.

TRUE = 1;
FALSE = 0;

% Constants used to switch between different test functions.
% Right now using EXP_0 to match the 0 Dirichlet BCs.
EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLYNOMIAL = 4;
EXPONENT_0 = 5;
global test_solution
test_solution = EXPONENT_0;

% Constants used to switch between different test domains.
% Right now only the square domain is defined in 3D.
CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;
RECTANGLE = 5;
DIAMOND_2 = 6;
SQUARE = 7;
global domain
domain = SQUARE;

% Description of spatial grid.  We divide both the x and y dimensions of
% the problem into the same number of gridpoints.  First, identify the min 
% and max values in both directions.  The actual domain of the PDE will be 
% a subset of the big rectangle defined by these values.
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
end

% At this point we override the above min/max definitions to
% implement the method on a 3D unit cube.  The above definitions
% can be corrected later.

x_min = 0;
x_max = 1;
y_min = 0;
y_max = 1;
z_min = 0;
z_max = 1;


N = 80;
hx = (x_max - x_min)/N;
hy = (y_max - y_min)/N;
hz = (z_max - z_min)/N;

% For the simple problem, there is no need for a complicated data 
% structure.  The entire problem is defined on a rectangular prism
% with the boundary points of the domain on the edges of prism.
grid = zeros(N+1, N+1, N+1);
% Fill in the initial value of the PDE
for i = 0:N
   for j = 0:N
      for k = 0:N
         grid(i+1, j+1, k+1) = g1(h*i, h*j, h*k);
      end
   end
end

% Define temporal grid.  Use a timestep of the same size as the spatial 
% grids. Go to time 1.
tau = h;
M = 1;

% Allocate space for plots.
plot_data = zeros(N+1,N+1,N+1,9);

% Save initial state to plot for debugging
for i = 0:N 
    for j = 0:N
        for k = 0:N
            plot_data(i+1, j+1, k+1, 1) = grid(i+1, j+1, k+1);
        end
    end
end

% Grids used to hold partial steps for the method
V1 = zeros(N+1, N+1, N+1);
V2 = zeros(N+1, N+1, N+1);
V3 = zeros(N+1, N+1, N+1);

% Grid used for forcing terrm
load_array = zeros(N+1,N+1,N+1);

% Evaluate each timestep.
for m = 1:M
    % Populate the array that describes the effect of the forcing term
    for i = 1:N-1
       for j = 1:N-1
           for k = 1:N-1
            load_array(i,j) = (1/2)*(f(h*i, h*j, h*k, tau*(m-1)) + ...
                f(h*i, h*j, h*k, tau*m)); 
           end
       end
    end
    
    % Create an array of values for V1 at each interior gridpoint
    for i = 1:N-1
       for j = 1:N-1
          V(i,j) = U(i+1,j+1) + ...
              (tau / (2*h^2)) * (U(i+1,j) - 2*U(i+1,j+1) + U(i+1,j+2));
       end
    end
    
    
    % Calculate the value of the function V1 at every interior gridpoint.
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
                (tau / (col(c).btf.h_prime + hy) ) * ...
                ( (grid(i,j+1).U - grid(i,j).U) / hy ...
                - (grid(i,j).U - col(c).btf.U) / ...
                col(c).btf.h_prime );
            j = j + 1;
            while (j < col(c).j_max)
                grid(i,j).V = grid(i,j).U + ...
                    (tau / (2*hy^2)) * (grid(i,j+1).U ...
                    - 2 * grid(i,j).U ...
                    + grid(i,j-1).U);
                j = j + 1;
            end
            grid(i,j).V = grid(i,j).U + ...
                (tau / (hy + col(c).btl.h_prime) ) * ...
                ( (col(c).btl.U - grid(i,j).U) / ...
                col(c).btl.h_prime ...
                - (grid(i,j).U - grid(i,j-1).U) / hy );                
       end
    end
    
    % Plot V for debuggging
    if m == 1
       for i = 1:N - 1
          for j = 1:N - 1
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
            b(1) = b(1) + (tau/((hx + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((hx + btl.h_prime) * btl.h_prime) ) * btl.U;
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
       x = A\b;
       
       % Store the values of U for the current row
       for i = i_min:i_max
           grid(i,j).U = x(i - i_min + 1);
       end                   
    end
    
    
    % Plot half-timestep solution for debugging
    %{
    if m == 1
       for i = 1:N - 1
          for j = 1:N - 1
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
    %}
    
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
            b(1) = b(1) + (tau/((hy + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((hy + btl.h_prime) * btl.h_prime) ) * btl.U;
       end
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*hy^2)) * A;
       % Adjust first and last row since boundary points are unevenly
       % spaced.  First consider case when A is 1 x 1
       if length == 1
            A(1,1) = 1 + tau / (btf.h_prime * btl.h_prime);
       % Then case when A is 2 x 2 or larger
       else
            A(1,1) = tau / (hy * btf.h_prime);
            A(1,2) = - tau / (hy * (hy + btf.h_prime));
            A(length, length - 1) = - tau / (hy * (hy + btl.h_prime));
            A(length, length) = tau / (hy * btl.h_prime);       
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
       for i = 1:N - 1
          for j = 1:N - 1
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
        for i = 1:N - 1
          for j = 1:N - 1
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

subplot(3,3,4)
surf(X,Y,plot_data(:,:,4))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after one timestep')

subplot(3,3,5)
surf(X,Y,plot_data(:,:,5))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after one timestep')

subplot(3,3,6)
surf(X,Y,plot_data(:,:,6));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after one timestep')

subplot(3,3,7)
surf(X,Y,plot_data(:,:,7))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after final timestep')

subplot(3,3,8)
surf(X,Y,plot_data(:,:,8))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after final timestep')

subplot(3,3,9)
surf(X,Y,plot_data(:,:,9));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after final timestep')

end