function [polar_coordinates] = custom_cart2pol(A)
    col_norm = sqrt(sum(A.^2,1));
    col_bearing = atan2(A(2,:), A(1,:));
    polar_coordinates = [col_norm; col_bearing];
end

