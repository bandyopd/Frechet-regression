function [dcov] = compute_dcov(X, Y)
% compute_dcov compute the dcov between X and Y
%  X: n*p matrix, each row is a sample
%  Y: n*3 matrix, each row is a sample
[n, p] = size(X);
dcov = zeros(p, 1);
for i = 1:p
    x_i = X(:, i);
    A = zeros(n, n);
    B = zeros(n, n);
    for j = 1:n
        for k = 1:n
            A(j, k) = abs(x_i(j) - x_i(k));
            tmp = sum(Y(j, :).*Y(k, :));
            if abs(tmp) > 1
                B(j, k) = acos(round(tmp));
            else
                B(j, k) = acos(tmp);
            end
        end
    end
    % row mean and col mean
    A_row_mean = mean(A, 2);
    A_col_mean = mean(A, 1); 
    B_row_mean = mean(B, 2);
    B_col_mean = mean(B, 1);
    A_mean = mean(A, "all");
    B_mean = mean(B, "all");
    for j = 1:n
        for k = 1:n
            A(j, k) = A(j, k) - A_row_mean(j) - A_col_mean(k) + A_mean;
            B(j, k) = B(j, k) - B_row_mean(j) - B_col_mean(k) + B_mean;
        end
    end
    dcov(i) = sum(A.*B, "all") / n^2;
end

end