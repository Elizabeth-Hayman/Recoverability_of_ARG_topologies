function p_matrix_sol = solve_full_tree_topology(n, N0, theta)

% Solves the recursion for the probabilities that the full tree topology is
% known, starting from a sample size n, with a mutation rate theta, and
% allowing up to N0 skipped mutations.
% p(n0, nl, nf) is the probability conditioned on nl lineages remaining,
% with nf in state 2, with 

p_matrix = zeros(N0+1, n, n+1);
% Set up matrix for recursions, initialise as tree terminates at nl=1, 
% where it is trivially fully detectable for any N0.
p_matrix(:, 1, 1) = ones(N0+1, 1);
t=theta/2; %more convenient to work with

p_sol = zeros(N0+1, n); %solution matrix
p_sol(:,1) = ones(N0+1, 1);

%The tree prob reutrn function returns the probability if it exists (the
%indices make physical sense) and 0 if else.
for n0 = 0:N0
    for nl = 2:n
        for nf = fliplr(0:nl)
            p_matrix(n0+1, nl, nf+1) =(nf/2*(nf-1)*tree_prob_return(p_matrix, n0, nl-1, nf-2, n, N0)+...
                (nl-nf)*t*tree_prob_return(p_matrix, n0, nl, nf+1, n, N0)+...
                nf*(nl-nf)*tree_prob_return(p_matrix, n0-1, nl-1, nf-1, n, N0)+...
                (nl-nf)*(nl-nf-1)/2*tree_prob_return(p_matrix, n0-2, nl-1, nf, n, N0)) / (nl/2*(nl-1)+ t*(nl-nf));
        end
        p_sol(n0+1, nl) = p_matrix(n0+1, nl, nl);
    end
end

% Final probability given by p_matrix(N0, nl=n, nf=n) 
p_matrix_sol = p_sol;
end
