clf
tic
R=2;
N0=0;
n=50;
rho=0.1;
THETA=[10,10,1000,10^10];
z=.5;

for theta= THETA
    P = solve_full_topology(R,n,N0,rho,theta,z);
    [~, q] = solve_number_recombs(n, R, rho);
    for n0 = 0:N0
        probs = zeros(1, n);
        
        for nl=1:n
            if q(nl,1) ==0
                probs(nl) = 1;
            else
                probs(nl) = P(n0+1, nl, nl+1, 1, 1, 1, 1)/q(nl,1);
            end
        end
        plot(probs, 'LineWidth',3);
        hold on
        xlim([0,n])
        ylim([0,.2]);
    end
end
legend('\theta=10', '\theta=100', '\theta=1000', '\theta=\infty')
    %legend('n_0=0', 'n_0=1', 'n_0=2', 'n_0=3')
toc