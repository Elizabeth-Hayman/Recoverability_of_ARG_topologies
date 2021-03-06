THETA=[1, 10, 100];
n = 50;
N0=3;

%Plot probabilties that we have R recombs total for range of small R
for i = 1:length(THETA)
    theta= THETA(i);
    figure(i);clf;
    p_matrix = solve_full_tree_topology(n, N0, theta);
    
    for n0=0:N0
        plot(p_matrix(n0+1, :), 'LineWidth', 3);
        hold on
    end
    legend("N_0=0","N_0=1","N_0=2","N_0=3")
end