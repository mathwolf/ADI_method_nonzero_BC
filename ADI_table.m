function ADI_table()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions. This program tests the method using the function
% e^(2x+3y+t).  It uses a range of grid sizes and displays
% the associated error measurements in a table.

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the grid 
% size.

for p = 1:5
    % Use a spatial grid of 0.1 times 2 to the power p-1
    % Start from the assumption that our test function is defined on the
    % domain [-1,1] in both x and y
    h = 1./(100 + 50*(p-1));
    N = 2./h + 1;
    x_min = -1.;
    y_min = -1.;
    stencil = make_stencil(h, N, x_min, y_min);
    
    % Set up data structure that is defined only on gridpoints that are
    % interior to the domain of the problem
    for i = 1:N
        % Go along the columns starting at the left
        for j = 1:N
            if stencil(i,j) == 0
                % Not an interior point
                grid(i,j).on = 0; % FALSE
            else
                % Interior point
                grid(i,j).on = 1; % TRUE
                grid(i,j).x = x_min + h * (i-1);
                grid(i,j).y = y_min + h * (j-1);
                grid(i,j).U = g1(grid(i,j).x, grid(i,j).y); % using IC
                grid(i,j).V = 0.;
            end
        end
    end
    
    % Set up rows for the data structure
    n_rows = 0;
    for j = 1:N
        % Go along row until we reach the end or find an interior point
        i = 1;
        while (i <= N) && (grid(i,j).on == 0) % FALSE
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
        % Go along row until we reach the end or find an exterior point
        while (i <= N) && (grid(i,j).on == 1) % TRUE
            i = i + 1; 
        end
        row(n_rows).endi = i - 1;
    
        % Add boundary value terms to each nonempty row
        row(n_rows).btf.y = y_min + h * (j-1);
        row(n_rows).btf.x = psi1(row(n_rows).btf.y);
        row(n_rows).btf.h_prime = x_min + h * (row(n_rows).starti - 1) ...
            - row(n_rows).btf.x;
    
        row(n_rows).btl.y = y_min + h * (j-1);
        row(n_rows).btl.x = psi2(row(n_rows).btl.y);
        row(n_rows).btl.h_prime = row(n_rows).btl.x ...
            - ( x_min + h * (row(n_rows).endi - 1) );   
    end

    % Set up columns on data structure
    n_cols = 0;
    for i = 1:N
        % Go through col until we reach the end or find an interior point
        j = 1;
        while (j <= N) && (grid(i,j).on == 0) % FALSE
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
        % Go through col until we reach the end or the next exterior point
        while (j <= N) && (grid(i,j).on == 1) % TRUE
            j = j + 1; 
        end
        col(n_cols).endj = j - 1;
    
        % Add boundary value terms to each nonempty col
        col(n_cols).btf.x = x_min + h * (i-1);
        col(n_cols).btf.y = phi1(col(n_cols).btf.x);
        col(n_cols).btf.h_prime = y_min + h * (col(n_cols).startj - 1) ...
            - col(n_cols).btf.y;    
        col(n_cols).btf.U = ... % use initial condition
            g1(col(n_cols).btf.x, col(n_cols).btf.y);

        col(n_cols).btl.x = x_min + h * (i-1);
        col(n_cols).btl.y = phi2(col(n_cols).btl.x);
        col(n_cols).btl.h_prime = col(n_cols).btl.y ...
            - ( y_min + h * (col(n_cols).endj - 1) );
        col(n_cols).btl.U = ... % use initial condition
            g1(col(n_cols).btl.x, col(n_cols).btl.y);
        
    end

    % Use a temporal gridsize of the same size as the spatial size, go to 
    % 1 second
    tau = h;
    M = 1./tau;

    % Loop to go through timesteps
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
            btf.U = g2(row(r).btf.x, row(r).btf.y, tau*(m-0.5));
            btf.h_prime = row(r).btf.h_prime;
            btl.U = g2(row(r).btl.x, row(r).btl.y, tau*(m-0.5));
            btl.h_prime = row(r).btl.h_prime;       
       
            % Add effect of boundary pts to rhs of matrix equation
            b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       
            % Set up tridiagonal matrix on lhs of equation
            A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
            A = (tau/(2*h^2)) * A;
            % Adjust first and last row since boundary points are unevenly
            % spaced
            A(1,1) = tau / (h * btf.h_prime);
            A(1,2) = - tau / (h * (h + btf.h_prime));
            A(length, length - 1) = - tau / (h * (h + btl.h_prime));
            A(length, length) = tau / (h * btl.h_prime);
            A = speye(length) + A ;
       
            % Solve the vector equation
            x = A\b;
       
            % Update the approximate solution for the current row
            for i = starti:endi
                grid(i,j).U = x(i - starti + 1);
            end                   
        end

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
            btf.U = g2(col(c).btf.x, col(c).btf.y, tau*m);
            btf.h_prime = col(c).btf.h_prime;
            btl.U = g2(col(c).btl.x, col(c).btl.y, tau*m);
            btl.h_prime = col(c).btl.h_prime;
            
            col(c).btf.U = btf.U;
            col(c).btl.U = btl.U;
       
            % Add effect of boundary pts to rhs of equation
            b(1) = b(1) + (tau/((h + btf.h_prime) * btf.h_prime) ) * btf.U;
            b(length) = b(length) + ...
                (tau/((h + btl.h_prime) * btl.h_prime) ) * btl.U;
       
            % Set up tridiagonal matrix on lhs of equation
            A = spdiags([-ones(length,1), 2*ones(length,1), -ones(length,1)], ...
                        [-1,0,1], length, length);
            A = (tau/(2*h^2)) * A;
            % Adjust first and last row since boundary points are unevenly
            % spaced
            A(1,1) = tau / (h * btf.h_prime);
            A(1,2) = - tau / (h * (h + btf.h_prime));
            A(length, length - 1) = - tau / (h * (h + btl.h_prime));
            A(length, length) = tau / (h * btl.h_prime);       
            A = speye(length) + A ;
       
            % Solve the vector equation
            x = A\b;

            % Update the approximate solution for the current row
            for j = startj:endj
                grid(i,j).U = x(j - startj + 1);
            end
        end
    end
    
    % Evaluate error at final timestep
    max_nodal_error = 0;
    nodal_error_squared = 0;    
    root_mean_sq_error = 0;

    for r = 1:n_rows
        j = row(r).j;
        starti = row(r).starti;
        endi = row(r).endi;
        for i = starti:endi
           current_error = abs(grid(i,j).U - u(grid(i,j).x, grid(i,j).y, 1));
           if current_error > max_nodal_error
              max_nodal_error = current_error; 
           end
           nodal_error_squared = nodal_error_squared + current_error^2;
        end
    end
    root_mean_sq_error = sqrt(h^2 * nodal_error_squared);

    
    % Write data to table
    table_data(p,1) = h;
    table_data(p,2) = max_nodal_error;
    if p ~= 1
       % estimate order of convergence for max err
       table_data(p,3) = log(table_data(p-1,2) / table_data(p,2)) / ...
           log(table_data(p-1,1) / table_data(p,1)); 
    end   
    table_data(p,4) = h;
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
