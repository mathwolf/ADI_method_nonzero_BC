function ADI_graph()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N= 9;
h = 0.25;
x_min = -1.;
y_min = -1.;

stencil = make_stencil(h, N, x_min, y_min);

tau = h;
M = 1./tau;


% Set up matrix of gridpoints to represent domain of problem
for i = 1:N
   % Go along the columns starting at the left
   for j = 1:N
       if stencil(i,j) == 0
           % Not an interior point
           interior_pts(i,j).on = 0; % FALSE
       else
           % Interior point
           interior_pts(i,j).on = 1; % TRUE
           interior_pts(i,j).x = x_min + h * (i-1);
           interior_pts(i,j).y = y_min + h * (j-1);
           interior_pts(i,j).U = ... % use initial condition
               g1(interior_pts(i,j).x, interior_pts(i,j).y);
           interior_pts(i,j).V = 0.;
       end
   end
end

% Set up row structure for matrix equations
n_rows = 0;
for j = 1:N
    % Go along row until we reach the end or find an interior point
    i = 1;
    while (i <= N) && (interior_pts(i,j).on == 0) % FALSE
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
    while (i <= N) && (interior_pts(i,j).on == 1) % TRUE
       i = i + 1; 
    end
    row(n_rows).endi = i - 1;
    
    % Add boundary value terms to each nonempty row
    row(n_rows).btf.y = y_min + h * (j-1);
    row(n_rows).btf.x = psi1(row(n_rows).btf.y);
    row(n_rows).btf.h_prime = x_min + h * (row(n_rows).starti - 1) ...
        - row(n_rows).btf.x;
    row(n_rows).btf.U = ... % use initial condition
        g1(row(n_rows).btf.x, row(n_rows).btf.y);
    row(n_rows).btl.y = y_min + h * (j-1);
    row(n_rows).btl.x = psi2(row(n_rows).btl.y);
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
    while (j <= N) && (interior_pts(i,j).on == 0) % FALSE
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
    while (j <= N) && (interior_pts(i,j).on == 1) % TRUE
       j = j + 1; 
    end
    col(n_cols).endj = j - 1;
    
    % Add boundary value terms to each nonempty col
    % Add boundary value terms to each nonempty row
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

% Now that the spatial domain has been established, 
% evaluate each timestep 
for m = 1:2
    
    % Calculate the value of the function V at every interior gridpoint
    for c = 1:n_cols
       i = col(c).i;
       j = col(c).startj;
       
       % First, consider special case where our column holds only 1
       % interior point
       if j == col(c).endj
           interior_pts(i,j).V = ... % formula using (5)
                interior_pts(i,j).U + (tau / (2*h)) * ...
                ( (col(c).btl.U - interior_pts(i,j).U) / ...
                    col(c).btl.h_prime ...
                    - (interior_pts(i,j).U - col(c).btf.U) / ...
                    col(c).btf.h_prime );
           
       % Otherwise the row holds two or more interior points
       else
            interior_pts(i,j).V = ... % formula for first point, using (5)
                interior_pts(i,j).U + (tau / (2*h)) * ...
                ( (interior_pts(i,j+1).U - interior_pts(i,j).U) / h ...
                    - (interior_pts(i,j).U - col(c).btf.U) / ...
                    col(c).btf.h_prime );
            j = j + 1;
            while (j < col(c).endj)
                interior_pts(i,j).V = interior_pts(i,j).U + ...
                    (tau / (2*h^2)) * (interior_pts(i,j+1).U ...
                        - 2 * interior_pts(i,j).U ...
                        + interior_pts(i,j-1).U);
                j = j + 1;
            end
            interior_pts(i,j).V = ... % formula for last point, using (5)
                interior_pts(i,j).U + (tau / (2*h)) * ...
                ( (col(c).btl.U - interior_pts(i,j).U) / ...
                    col(c).btl.h_prime ...
                - (interior_pts(i,j).U - interior_pts(i,j-1).U) / h );
                
            
                ( (col(c).btf.U - interior_pts(i,j).U) / ...
                    col(c).btf.h_prime ...
                 - (interior_pts(i,j).U - interior_pts(i,j+1).U) / h );
                
       end
    end
    
end

end
