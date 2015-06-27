function ADI_nonzero_BC_table()

% Matlab R2013a

% ADI method for solution of parabolic PDE with Dirichlet boundary
% conditions. This program tests the method using the function
% e^(x+y+t).  It uses a range of grid sizes and displays
% the associated error measurements in a table.

% Table for storing error data
table_data = zeros(5,6);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.1 times 2 to the power p-1
    N = 10 * 2^(p-1);
    h = 1./N;
    
    % Use a temporal spacing of the same size, go to 1 second
    M = N;
    tau = 1./M;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = zeros(N+1,N+1);
    for i = 0:N
        for j = 0:N
            U(i+1,j+1) = g1(h*i, h*j); 
        end
    end

    % Tridiagonal matrices used on LHS of ADI method
    tridiagonal = spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
    A = speye(N-1) + (tau/(2*h^2))*tridiagonal;

    % Matrix V used on RHS of matrix formulas, interior points only
    V = zeros(N-1,N-1);

    % Array for forcing term
    load_vector = zeros(N-1,N-1);

    % Boundary values used on RHS of matrix formulas
    boundary_term = zeros(N-1,1);

    % Step through the scheme
    for m = 1:M
        % Create an array for the value of the forcing term at each 
        % gridpoint at one half of a timestep
        for i = 1:N-1
            for j = 1:N-1
                load_array(i,j) = (1/2)*(f(h*i, h*j, tau*(m-1)) + ...
                    f(h*i, h*j, tau*m)); 
            end
        end

        % Create an array of values for V at each interior gridpoint
        for i = 1:N-1
            for j = 1:N-1
                V(i,j) = U(i+1,j+1) + (tau / (2*h^2)) * ...
                    (U(i+1,j) - 2*U(i+1,j+1) + U(i+1,j+2));
            end
        end
        
        % Find the approximte vlues at the half timestep
        % Need to solve system once for each j = 2,...,N
        for j = 1:N-1  
            % Calculate the boundary values for the first and last points
            btf = (1/2) * (g2(0, h*j, tau*m) + g2(0, h*j, tau*(m-1))) ...
                - (tau/(4*h^2)) * (g2(0, h*(j-1), tau*m) - ...
                2*g2(0, h*j, tau*m) + g2(0, h*(j+1), tau*m)) ...
                + (tau/(4*h^2)) * (g2(0, h*(j-1), tau*(m-1)) - ...
                2*g2(0, h*j, tau*(m-1)) + g2(0, h*(j+1), tau*(m-1)));
            btl = (1/2) * (g2(1, h*j, tau*m) + g2(1, h*j, tau*(m-1))) ...
                - (tau/(4*h^2)) * (g2(1, h*(j-1), tau*m) - ...
                2*g2(1, h*j, tau*m) + g2(1, h*(j+1), tau*m)) ...
                + (tau/(4*h^2)) * (g2(1, h*(j-1), tau*(m-1)) - ...
                2*g2(1, h*j, tau*(m-1)) + g2(1, h*(j+1), tau*(m-1)));
            % Find rhs of vector equation
            b = V(1:N-1,j) + (tau/2) * load_array(1:N-1,j);
            % Add effect of boundary values of rhs
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve vector equation to get points for half timestep
            U(2:N,j+1) = A\b;
        end
    
        % Now find the interior points for next whole timestep
        % Solve system one for each i = 2, ... , N
        for i = 1:N-1 
            % Find boundary values for first and last points
            btf = g2(h*i, 0, tau*m);
            btl = g2(h*i, h*N, tau*m);
            % Calculate rhs of next vector equation
            b = 2*U(i+1,2:N)' - V(i,1:N-1)';
            % Add effect of boundary values
            b(1) = b(1) + (tau/(2*h^2)) * btf;
            b(N-1) = b(N-1) + (tau/(2*h^2)) * btl;
            % Solve vector equation to get points for next full timestep
            U(i+1,2:N) = (A\b)';       
        end
    
        % Update boundary points
        for i = 0:N
            U(1,i+1) = g2(0, h*i, tau*m); 
            U(N+1,i+1) = g2(h*N, h*i, tau*m);
            U(i+1,1) = g2(h*i, 0, tau*m);
            U(i+1,N+1) = g2(h*i, h*N, tau*m);
        end

    end
    
    % Evaluate error at interior gridpoints at final timestep
    max_nodal_error = 0;
    nodal_error_squared = 0;    
    root_mean_sq_error = 0;

    for i = 1:N-1
       for j = 1:N-1
           current_error = abs(U(i+1,j+1) - u(h*i, h*j, 1));
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

