function ADI_nonzero_BC_table()
% Summary of this function goes here
%   Detailed explanation goes here

% Table for storing error data
table_data = zeros(5,5);

% Check five different grid sizes. Each step will decrease the size by 
% half.

for p = 1:5
    
    % Use a spatial grid of 0.1 times 2 to the power p-1
    N = 10 * 2^(p-1) + 1;
    h = 1./(N-1);
    
    % Use a temporal spacing of the same size, go to 1 second
    tau = h;
    M = floor(1./tau) + 1;

    % Permittivity parameter for material
    a = 1.;

    % Set up matrix of gridpoints
    % Use g1 for initial conditions
    U = sparse(N,N);
    for i = 1:N
        for j = 1:N
            U(i,j) = g1(h*(i-1), h*(j-1)); 
        end
    end



    % Matrices used on LHS of ADI method
    Bh = (1/h^2) * spdiags([-ones(N-2,1), 2*ones(N-2,1), -ones(N-2,1)], ...
                        [-1,0,1], N-2, N-2);
    R = chol(speye(N-2) + (a*tau/2)*Bh);

    % Matrix used to store alternating step
    U_twiddles = sparse(N,N);

    % Array for forcing term
    load_vector = sparse(N-2,N-2);

    boundary_term = sparse(N-2,1);

    % Step through the scheme
    for m = 1:M
        % Create an array for the value of the forcing term at each 
        % gridpoint at one half of a timestep
        for i = 1:N
            for j = 1:N
                load_vector(i,j) = (1/2)* ...
                    (f(h*(i-1), h*(j-1), tau*(m-1)) + ...
                    f(h*(i-1), h*(j-1), tau*m)); 
            end
        end
        
        % Find the boundary values for the alternating direction.  Only the
        % values along the x boundary are needed by the formulas.
        for j = 2:N-1
            U_twiddles(1,j) = -a*tau/(4*h^2) * g2(0, h*(j-2), tau*m) + ...
                (1/2 + a*tau/(2*h^2)) * g2(0, h*(j-1), tau*m) + ...
                -a*tau/(4*h^2) * g2(0, h*j, tau*m) + ...
                a*tau/(4*h^2) * g2(0, h*(j-2), tau*(m-1)) + ...
                (1/2 - a*tau/(2*h^2)) * g2(0, h*(j-1), tau*(m-1)) + ...
                a*tau/(4*h^2) * g2(0, h*j, tau*(m-1));
            U_twiddles(N,j) = -a*tau/(4*h^2) * g2(1, h*(j-2), tau*m) + ...
                (1/2 + a*tau/(2*h^2)) * g2(1, h*(j-1), tau*m) + ...
                -a*tau/(4*h^2) * g2(1, h*j, tau*m) + ...
                a*tau/(4*h^2) * g2(1, h*(j-2), tau*(m-1)) + ...
                (1/2 - a*tau/(2*h^2)) * g2(1, h*(j-1), tau*(m-1)) + ...
                a*tau/(4*h^2) * g2(1, h*j, tau*(m-1));       
        end
    

    
        % Find the alternate matrix U_twiddles 
        % Need to solve system once for each i = 2,...,N-1
        for i = 2:N-1  
        
            boundary_term(1,1) = U_twiddles(1,i);
            boundary_term(N-2,1) = U_twiddles(N,i);
            b = (a*tau)/(2*h^2) * U(2:N-1,i-1) + ...
                (1 - 2*(a*tau)/(2*h^2)) * U(2:N-1,i) + ...
                (a*tau)/(2*h^2) * U(2:N-1,i+1) + ...
                (tau/2)*load_vector(2:N-1,i) + ...
                (a*tau)/(2*h^2) * boundary_term(:,1);

        %{
        b = (speye(N-2) - (a*tau/2)*Bh) * U(i,:)' + ...
            (tau/2)*load_vector(i,:)';
        %}
        %{
        if m == 1
           disp(full((speye(N-2) - (a*tau/2)*Bh) * U(i,:)'));
           disp(full(b)); 
        end
            %}
            b = R'\b;
        %{
        if m == 1
           disp(full(b)); 
        end
        %}
            U_twiddles(2:N-1,i) = R\b;
        %{
        if m ==1
           disp(full(U_twiddles(:,i))); 
        end
        %}
        end
        
        % Now find the interior points for next U
        % Solve system one for each j = 2, ... , N-1
        for j = 2:N-1
        
            boundary_term(1,1) = g2(h*(j-1), 0, tau*m) - U(j,1);
            boundary_term(N-2,1) = g2(h*(j-1), 1, tau*m) - U(j,N);
            b = 2*U_twiddles(j,2:N-1)' - ...
                (speye(N-2) - (a*tau/2)*Bh) *  U(j,2:N-1)' + ...
                (a*tau)/(2*h^2) * boundary_term(:,1);
       %{
       if m == 1
          disp(full(b));
       end
       %}
            b = R'\b;
            U(j,2:N-1) = (R\b)';       
       %{
       if m == 1
          disp(full(U(i,:))); 
       end
       %}
        end
    
        % Update boundary points
        for i = 1:N
            U(1,i) = g2(0, h*(i-1), tau*m); 
            U(N,i) = g2(1, h*(i-1), tau*m);
            U(i,1) = g2(h*(i-1), 0, tau*m);
            U(i,N) = g2(h*(i-1), 1, tau*m);
        end
        
    end
    
    % Evaluate error at final timestep
    max_nodal_error = 0;
    mean_sq_error = 0;
    nodal_error_squared = 0;
    
    for i = 2:N-1
       for j = 2:N-1
           current_error = abs(U(i,j) - u(h*(i-1), h*(j-1), tau*M));
           if current_error > max_nodal_error
              max_nodal_error = current_error; 
           end
           nodal_error_squared = nodal_error_squared + current_error^2;
       end
    end
    
    mean_sq_error = sqrt(nodal_error_squared / (N-2)^2);
    
    % Write data to table
    table_data(p,1) = h;
    table_data(p,2) = max_nodal_error;
    if p ~= 1
       table_data(p,3) = log2(table_data(p-1,2) / table_data(p,2)); 
    end
    
    table_data(p,4) = mean_sq_error;
    if p ~= 1
       table_data(p,5) = log2(table_data(p-1,4) / table_data(p,4)); 
    end
end
disp('    h               tau                     max error                       rms errror                            order of conv');
disp(table_data);

end

