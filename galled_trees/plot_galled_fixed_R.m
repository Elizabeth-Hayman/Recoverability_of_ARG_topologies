%Plot prob tree is galled conditional on a fixed number of recombinations
%in the history. Uniform prior taken over rho.
figure(1);clf;
n = 200;
m=20; %nmuber of data points across the range for rho
R_range = 1:5;
RHO=linspace(0, 1, 20);

for R = R_range
    cumsum_p = zeros(1,n);
    cumsum_q = zeros(1,n);
    for rho = RHO
        [~, p_matrix] = solve_number_galled(n, R, rho);
        cumsum_p = cumsum_p + p_matrix(1:n, 1, 1);
        [~, q_matrix] = solve_number_recombs(n, R, rho);
        cumsum_q = cumsum_q + q_matrix(1:n, 1);   
    end
    cumsum_p = cumsum_p/m;
    cumsum_q = cumsum_q/m;
    plot(cumsum_p(1:n,1)./cumsum_q(1:n,1), 'LineWidth', 3);
    hold on
end
ylim([0,1])
legend("R=1", "R=2", "R=3", "R=4", "R=5");
