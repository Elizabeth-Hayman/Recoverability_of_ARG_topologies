R_max= 3;
R_tot = 0:R_max;
rho = 0.1;
n = 250;
conditional = zeros(R_max+1, n);

%Plot probabilties that we have R recombs total for range of small R
figure(1);clf;
cumsum = zeros(1,n);
for R = R_tot
    [~, q_mat] = solve_number_recombs(n, R, rho);
    cumsum = cumsum + q_mat(1:n, 1);
    conditional(R+1, :) = q_mat(1:n, 1);
    plot(q_mat(1:n, 1), 'LineWidth', 3);
    hold on
end
plot(cumsum(1:n,1), 'LineWidth', 3);
hold on
legend("R=0","R=1","cumsum")
title("Total number of recombinations")

%Plot we have R galled recombs
%Plot probabilties for range of small R
figure(2);clf;
cumsum = zeros(1,n);
for R = R_tot
    [~, p_matrix] = solve_number_galled(n, R, rho);
    cumsum = cumsum + p_matrix(1:n, 1, 1);
    conditional(R+1, :) = conditional(R+1, :) ./ p_matrix(1:n, 1, 1);
    plot(p_matrix(1:n, 1, 1), 'LineWidth', 3);
    hold on
end
plot(cumsum(1:n,1), 'LineWidth', 3);
hold on
legend("R=0","R=1","cumsum");
title("Total number of galled recombinations");


%Plot we have a galled tree = sum_R(galled)/ sum_R(not galled)
%Plot probabilties for range of small R
figure(3);clf;
for R = R_tot
    plot(p_matrix(1:n, 1), 'LineWidth', 3);
    hold on
end
plot(cumsum(1:n,1), 'LineWidth', 3);
hold on
legend("R=0","R=1","cumsum");
title("Total number of galled recombinations");
