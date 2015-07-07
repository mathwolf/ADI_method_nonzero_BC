function ADI_graph()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

TRUE = 1;
FALSE = 0;

EXPONENT_1 = 1;
EXPONENT_2 = 2;
TRIG = 3;
POLY = 4;

global solution
solution = POLY;

CIRCLE = 1;
ELLIPSE = 2;
DIAMOND = 3;
ELL = 4;

domain = DIAMOND;

N= 21;
h = 0.1;
x_min = -1.;
y_min = -1.;

% Set up matrix of gridpoints to represent domain of problem
for i = 1:N
   % Go along the columns starting at the left
   for j = 1:N
       x = x_min + h * (i-1);
       y = y_min + h * (j-1);
       if (phi1(x, domain) < y) && (y < phi2(x, domain)) && ...
                (psi1(y, domain) < x) && (x < psi2(y, domain))
           % Interior point
           grid(i,j).on = TRUE; 
           grid(i,j).x = x;
           grid(i,j).y = y;
           grid(i,j).U = g1(x,y);  % use initial condition               
           grid(i,j).V = 0.;
       else
           % Exterior point
           grid(i,j).on = FALSE;
       end
   end
end

% Set up row structure for matrix equations
n_rows = 0;
for j = 1:N
    % Go along row until we reach the end or find an interior point
    i = 1;
    while (i <= N) && (grid(i,j).on == FALSE)
       i = i + 1; 
    end

    if i > N
       % We are at the end of an empty row, go to the next row
       continue;
    end
    
    % We are at the first interior point in this row
    n_rows = n_rows + 1;
    row(n_rows).j = j;
    row(n_rows).starti = i;
    while (i <= N) && (grid(i,j).on == TRUE)
       i = i + 1; 
    end
    row(n_rows).endi = i - 1;
    
    % Add boundary value terms to each nonempty row
    row(n_rows).btf.y = y_min + h * (j-1);
    row(n_rows).btf.x = psi1(row(n_rows).btf.y, domain);
    row(n_rows).btf.h_prime = x_min + h * (row(n_rows).starti - 1) ...
        - row(n_rows).btf.x;
    row(n_rows).btf.U = ... % use initial condition
        g1(row(n_rows).btf.x, row(n_rows).btf.y);
    
    row(n_rows).btl.y = y_min + h * (j-1);
    row(n_rows).btl.x = psi2(row(n_rows).btl.y, domain);
    row(n_rows).btl.h_prime = row(n_rows).btl.x ...
        - ( x_min + h * (row(n_rows).endi - 1) );   
    row(n_rows).btl.U = ... % use initial condition
        g1(row(n_rows).btl.x, row(n_rows).btl.y);
end

% Set up col structure for matrix equations
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
    
    % We are at the first interior point in this col
    n_cols = n_cols + 1;
    col(n_cols).i = i;
    col(n_cols).startj = j;
    while (j <= N) && (grid(i,j).on == TRUE) 
       j = j + 1; 
    end
    col(n_cols).endj = j - 1;
    
    % Add boundary value terms to each nonempty col
    col(n_cols).btf.x = x_min + h * (i-1);
    col(n_cols).btf.y = phi1(col(n_cols).btf.x, domain);
    col(n_cols).btf.h_prime = y_min + h * (col(n_cols).startj - 1) ...
        - col(n_cols).btf.y;    
    col(n_cols).btf.U = ... % use initial condition
        g1(col(n_cols).btf.x, col(n_cols).btf.y);
    
    col(n_cols).btl.x = x_min + h * (i-1);
    col(n_cols).btl.y = phi2(col(n_cols).btl.x, domain);
    col(n_cols).btl.h_prime = col(n_cols).btl.y ...
        - ( y_min + h * (col(n_cols).endj - 1) );
    col(n_cols).btl.U = ... % use initial condition
        g1(col(n_cols).btl.x, col(n_cols).btl.y);
end

% Define timestep

tau = h;
M = 1./tau;


%{
disp(n_rows);
disp(n_cols);

for i = 1:N
    for j = 1:N
        if interior_pts(i,j).on == 1
            disp(interior_pts(i,j));
        end
    end
end


for j=1:n_rows
   disp(row(j)); 
   disp(row(j).btf);
   disp(row(j).btl);
end

for i=1:n_cols
   disp(col(i)); 
   disp(col(i).btf);
   disp(col(i).btl);
end
%}

plot_data = zeros(N,N,6);
% Save initial state for plotting
%{

for i = 1:N
    for j = 1:N
        if interior_pts(i,j).on == 1
            plot_data(i,j,1) = interior_pts(i,j).U; 
            plot_data(i,j,2) = u(interior_pts(i,j).x, ...
                interior_pts(i,j).y, 0);
            plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
        end
    end
end
%}

% Now that the spatial domain has been established, 
% evaluate each timestep 
for m = 1:M
    
    % Calculate the value of the function V at every interior gridpoint
    for c = 1:n_cols
       i = col(c).i;
       j = col(c).startj;
       
       % First, consider special case where the column holds only 1
       % interior point
       if j == col(c).endj
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
            while (j < col(c).endj)
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
  
    %{
    % Plot V for debuggging
    if m == 1
       for i = 1:N
          for j = 1:N
              if interior_pts(i,j).on == 1
                 plot_data(i,j,3) = interior_pts(i,j).V; 
              end
          end
       end
    end
    %}
    
    % Find the approximation at the half-timestep
    % Here we solve a matrix-vector equation once for each row
    for r = 1:n_rows
       j = row(r).j;
       starti = row(r).starti;
       endi = row(r).endi;
       length = endi - starti + 1;
       
       % Create vectors for rhs of equation
       V = zeros(length,1);
       load_vector = zeros(length,1);
       for i = starti:endi
          V(i - starti + 1) = grid(i,j).V; 
          load_vector(i - starti + 1) = ...
              f(grid(i,j).x, grid(i,j).y, tau*(m-0.5));
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
       b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
       b(length) = b(length) + ...
           (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*h^2)) * A;
       
       % Adjust first and last row since boundary points are unevenly
       % spaced.  First consider case when A is 1 x 1
       if length == 1
           A(1,1) = tau / (btf.h_prime * btl.h_prime);           
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
       
       % Update the approximate solution for the current row
       for i = starti:endi
           grid(i,j).U = x(i - starti + 1);
       end                   
    end
    
    %{
    % Store data for plotting on selected steps
    
    if m == 1
       for i = 1:N
          for j = 1:N
              if interior_pts(i,j).on == 1
                 plot_data(i,j,4) = interior_pts(i,j).U; 
                 plot_data(i,j,5) = u(interior_pts(i,j).x, ...
                     interior_pts(i,j).y, tau*(m-0.5));
                 plot_data(i,j,6) = ...
                     abs(plot_data(i,j,4) - plot_data(i,j,5));
              end
          end
       end
    end
    %}
    
    % Find interior points for the next whole timestep
    % Here we solve a matrix/vector equation once for each column
    for c = 1:n_cols
        i = col(c).i;
        startj = col(c).startj;
        endj = col(c).endj;
        length = col(c).endj - col(c).startj + 1;
        
        % Create vectors for rhs of equation
        U = zeros(length,1);
        V = zeros(length,1);
        for j = startj:endj
            U(j - startj + 1) = grid(i,j).U;
            V(j - startj + 1) = grid(i,j).V;
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
       b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
       b(length) = b(length) + ...
           (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       
       % Set up tridiagonal matrix on lhs of equation
       A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
       A = (tau/(2*h^2)) * A;
       % Adjust first and last row since boundary points are unevenly
       % spaced.  First consider case when A is 1 x 1
       if length == 1
            A(1,1) = tau / (btf.h_prime * btl.h_prime);
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

       % Update the approximate solution for the current row
       for j = startj:endj
           grid(i,j).U = x(j - startj + 1);
       end
    end
    
    % Store data for plotting on selected steps
    
    if m == 1
       for i = 1:N
          for j = 1:N
              if grid(i,j).on == 1
                 plot_data(i,j,1) = grid(i,j).U; 
                 plot_data(i,j,2) = u(grid(i,j).x, ...
                     grid(i,j).y, tau*m);
                 plot_data(i,j,3) = ...
                     abs(plot_data(i,j,1) - plot_data(i,j,2));
              end
          end
       end
    elseif m == M
        for i = 1:N
          for j = 1:N
              if grid(i,j).on == 1
                 plot_data(i,j,4) = grid(i,j).U; 
                 plot_data(i,j,5) = u(grid(i,j).x, ...
                     grid(i,j).y, tau*m);
                 plot_data(i,j,6) = ...
                     abs(plot_data(i,j,4) - plot_data(i,j,5));
              end
          end
       end

    end
    
end

% Create plots
X = linspace(-1., 1., N);
Y = linspace(-1., 1., N);

%{ 
%Extra plots for debuging
figure
subplot(3,3,1)
waterfall(X,Y,plot_data(:,:,1))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution at initial condition')

subplot(3,3,2)
waterfall(X,Y,plot_data(:,:,2))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution at initial condition')

subplot(3,3,3)
waterfall(X,Y,plot_data(:,:,3));
colormap winter;
xlabel('x')
ylabel('y')
title('V at initial condition')
%}
subplot(2,3,1)
waterfall(X,Y,plot_data(:,:,1))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after one timestep')

subplot(2,3,2)
waterfall(X,Y,plot_data(:,:,2))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after one timestep')

subplot(2,3,3)
waterfall(X,Y,plot_data(:,:,3));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after one timestep')

subplot(2,3,4)
waterfall(X,Y,plot_data(:,:,4))
colormap winter
xlabel('x')
ylabel('y')
title('Approximate solution after final timestep')

subplot(2,3,5)
waterfall(X,Y,plot_data(:,:,5))
colormap winter
xlabel('x')
ylabel('y')
title('Exact solution after final timestep')

subplot(2,3,6)
waterfall(X,Y,plot_data(:,:,6));
colormap winter;
xlabel('x')
ylabel('y')
title('Error after final timestep')

end
