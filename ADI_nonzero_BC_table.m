function ADI_nonzero_BC_table()

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
    tau = h;
    M = floor(1./tau) + 1;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = zeros(N+1,N+1);
    for i = 0:N
        for j = 0:N
            U(i+1,j+1) = g1(h*i, h*j); 
        end
    end

    % Tridiagonal matrices used on LHS of ADI method
    Bh = (1/h^2) * spdiags([-ones(N-1,1), 2*ones(N-1,1), -ones(N-1,1)], ...
                        [-1,0,1], N-1, N-1);
    A = speye(N-1) + (tau/2)*Bh;

    % Matrix V used on left hand side of formula
    V = zeros(N+1,N+1);

    % Array for forcing term
    load_vector = zeros(N+1,N+1);

    boundary_term = zeros(N-1,1);

    % Step through the scheme
    for m = 1:M
        % Create an array for the value of the forcing term at each 
        % gridpoint at one half of a timestep
        for i = 0:N
            for j = 0:N
                load_array(i+1,j+1) = (1/2)*(f(h*i, h*j, tau*(m-1)) + ...
                    f(h*i, h*j, tau*m)); 
            end
        end

        % Create an array of values for V at each interior gridpoint
        % and at gridpoints along the x-boundary
        for i = 0:N
            for j = 1:N-1
                V(i+1,j+1) = U(i+1,j+1) + (tau / (2*h^2)) * ...
                    (U(i+1,j) - 2*U(i+1,j+1) + U(i+1,j+2));
            end
        end
    
        
        
        % Find the boundary values for the alternating direction.  Only the
        % values along the x boundary are needed by the formulas.
        for j = 1:N-1
            U(1,j+1) = (1/2) * g2(0, h*j, tau*m) ...
                - (tau/(4*h^2)) * (g2(0, h*(j-1), tau*m) - ...
                2*g2(0, h*j, tau*m) + g2(0, h*(j+1), tau*m)) ...
                + (1/2) * V(1,j+1);
            U(N+1,j+1) = (1/2) * g2(h*N, h*j, tau*m) ...
                - (tau/(4*h^2)) * (g2(h*N, h*(j-1), tau*m) - ...
                2*g2(h*N, h*j, tau*m) + g2(h*N, h*(j+1), tau*m)) ...
                + (1/2) * V(N+1,j+1);
        end
    
        % Find the alternate matrix at the half timestep
        % Need to solve system once for each j = 2,...,N
        for j = 1:N-1  
            boundary_term(1,1) = U(1,j+1);
            boundary_term(N-1,1) = U(N+1,j+1);
            b = V(2:N,j+1) + ...
                (tau/2) * load_array(2:N,j+1) + ...
                (tau/(2*h^2)) * boundary_term(:,1);
            U(2:N,j+1) = A\b;
        end
    
        % Now find the interior points for next whole timestep
        % Solve system one for each i = 2, ... , N
        for i = 1:N-1        
            boundary_term(1,1) = g2(h*i, 0, tau*m);
            boundary_term(N-1,1) = g2(h*i, h*N, tau*m);
            b = 2*U(i+1,2:N)' - V(i+1,2:N)' + ...
                (tau/(2*h^2)) * boundary_term(:,1);
            U(i+1,2:N) = (A\b)';       
        end
    
        % Update boundary points
        for i = 0:N
            U(1,i+1) = g2(0, h*i, tau*m); 
            U(N+1,i+1) = g2(h*N, h*i, tau*m);
            U(i+1,1) = g2(h*i, 0, tau*m);
            U(i+1,N+1) = g2(h*i, h*N, tau*m);
        end

        %{
        % For debugging, store gridpoints to plot on selected steps
        if m == 1
            for i = 1:N+1
                for j = 1:N+1
                    plot_data(i,j,1) = U(i,j);
                    plot_data(i,j,2) = u(h*(i-1), h*(j-1), tau*m);
                    plot_data(i,j,3) = plot_data(i,j,1) - plot_data(i,j,2);
                end
            end
        elseif m == M
            for i = 1:N+1
                for j = 1:N+1
                    plot_data(i,j,4) = U(i,j);
                    plot_data(i,j,5) = u(h*(i-1), h*(j-1), tau*m);
                    plot_data(i,j,6) = plot_data(i,j,4) - plot_data(i,j,5);
                end
            end 
        end
        %}

    end
    
    % Evaluate error at interior gridpoints at final timestep
    max_nodal_error = 0;
    root_mean_sq_error = 0;
    nodal_error_squared = 0;    
    for i = 1:N+1
       for j = 1:N+1
           current_error = abs(U(i,j) - u(h*(i-1), h*(j-1), tau*m));
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
       table_data(p,3) = log(table_data(p-1,2) / table_data(p,2)) / ...
           log(table_data(p-1,1) / table_data(p,1)); 
    end   
    table_data(p,4) = h;
    table_data(p,5) = root_mean_sq_error;
    if p ~= 1
       table_data(p,6) = log2(table_data(p-1,5) / table_data(p,5));
    end
    
    %{
    % Display plots for debugging
    % Create plots
    X = linspace(0., 1., N+1);
    Y = linspace(0., 1., N+1);

    figure(p)
    subplot(2,3,1)
    surf(X,Y,plot_data(:,:,1))
    colormap winter
    xlabel('x')
    ylabel('y')
    title('Approximate solution after one timestep')

    subplot(2,3,2)
    surf(X,Y,plot_data(:,:,2))
    colormap winter
    xlabel('x')
    ylabel('y')
    title('Exact solution after one timestep')

    subplot(2,3,3)
    surf(X,Y,plot_data(:,:,3));
    colormap winter;
    xlabel('x')
    ylabel('y')
    title('Error after one timestep')

    subplot(2,3,4)
    surf(X,Y,plot_data(:,:,4))
    colormap winter
    xlabel('x')
    ylabel('y')
    title('Approximate solution after final timestep')

    subplot(2,3,5)
    surf(X,Y,plot_data(:,:,5))
    colormap winter
    xlabel('x')
    ylabel('y')
    title('Exact solution after final timestep')

    subplot(2,3,6)
    surf(X,Y,plot_data(:,:,6));
    colormap winter;
    xlabel('x')
    ylabel('y')
    title('Error after final timestep')
    %}

end

% Display table of error information
disp('    h              max error                         order of conv');
disp(table_data(:,1:3));
disp('    h              rms error                         order of conv');
disp(table_data(:,4:6));
end

