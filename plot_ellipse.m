% Function to plot ellipses
function plot_ellipse(x, y, color)
    P = cov([x y]); % Covariance matrix
    mean_x = mean(x);
    mean_y = mean(y);
    [eigvec, eigval] = eig(P);
    
    % Get the index of the largest eigenvector
    [largest_eigvec_ind_c, ~] = find(eigval == max(max(eigval)));
    largest_eigvec = eigvec(:, largest_eigvec_ind_c);
    
    % Get the largest eigenvalue
    largest_eigval = max(max(eigval));
    
    % Get the smallest eigenvector and eigenvalue
    if(largest_eigvec_ind_c == 1)
        smallest_eigval = max(eigval(:,2));
        smallest_eigvec = eigvec(:,2);
    else
        smallest_eigval = max(eigval(:,1));
        smallest_eigvec = eigvec(:,1);
    end
    
    % Calculate the angle between the x-axis and the largest eigenvector
    angle = atan2(largest_eigvec(2), largest_eigvec(1));
    
    % This angle is between -pi and pi.
    % Let's shift it such that the angle is between 0 and 2pi
    if(angle < 0)
        angle = angle + 2*pi;
    end
    
    % Get the coordinates of the data mean
    chisquare_val = 2.4477; % 95% confidence interval
    theta_grid = linspace(0,2*pi);
    phi = angle;
    X0=mean_x;
    Y0=mean_y;
    a=chisquare_val*sqrt(largest_eigval);
    b=chisquare_val*sqrt(smallest_eigval);
    
    % the ellipse in x and y coordinates 
    ellipse_x_r  = a*cos( theta_grid );
    ellipse_y_r  = b*sin( theta_grid );
    
    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    %let's rotate the ellipse to some angle phi
    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
    
    % Draw the error ellipse
    plot(r_ellipse(:,1) + X0, r_ellipse(:,2) + Y0,'-', 'Color', color);
end
