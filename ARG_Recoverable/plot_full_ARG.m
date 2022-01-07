%Current configuration will plot the probability of full ARG detection conditional on 2 recombinations for a variety of theta.
%Can be configured 
clf
tic

R=2; %configurable params
N0=0;
n=50;
rho=0.1;
THETA=[10,100,1000,10^10]; %selection of theta. 10^10 ~ infinity
z=.5;

for theta= THETA
    P = solve_full_topology(R,n,N0,rho,theta,z);
    [~, q] = solve_number_recombs(n, R, rho); %probability a sample has R recombination, used for conditioning

    probs = zeros(1, n); %initialise an empty probabilty vector each for loop
    for nl=1:n
        if q(nl,1) ==0 
            probs(nl) = 0; %if loop prevent NaN entries: if p and q are zero, want the ratio to also be 0
        else
            probs(nl) = P(n0+1, nl, nl+1, 1, 1, 1, 1)/q(nl,1); %calculate the conditional probability for non-zero entries
        end
    end
    plot(probs, 'LineWidth', 3);
    hold on
end
%Configure the graph
legend('\theta=10', '\theta=100', '\theta=1000', '\theta=\infty')
xlim([1,n])
ylim([0,.2]);

toc %outputs run time
