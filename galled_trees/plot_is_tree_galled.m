%Plot prob tree is galled, ie 
%cumulative probabilties summed to large R, for range of small rho
figure(1);clf;
n = 200;
R_range = 0:75; % sum to infinity by taking R suitably large
RHO=[0, .1, .5, 1, 10]; %rho datapoints
for rho = RHO
    cumsum_p = zeros(1,n);
    cumsum_q = zeros(1,n);
    for R = R_range
        [~, p_matrix] = solve_number_galled(n, R, rho);
        cumsum_p = cumsum_p + p_matrix(1:n, 1, 1);
        [~, q_matrix] = solve_number_recombs(n, R, rho);
        cumsum_q = cumsum_q + q_matrix(1:n, 1);
        
    end
    plot(cumsum_p(1:n,1)./cumsum_q(1:n,1), 'LineWidth', 3);
    hold on
end
legend("\rho=0", "\rho=0.1", "\rho=0.5", "\rho=1", "\rho=10")
