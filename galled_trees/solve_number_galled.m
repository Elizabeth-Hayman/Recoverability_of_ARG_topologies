function [p_sol, p_matrix_sol] = solve_number_galled(n, R, rho)
% Solves the recursion for the probabilities that a total of R
% recombinations occur in a sample history size n, all of which are galled.
% Returns a single probability and also whole matrix calculated in
% recursions.
% p(nl, r0, r) is the probability conditioned on r recombs happened so far
% in history, r0 recombs currently open, with nl lineages remaining

p_matrix = zeros(n+R, R+1, R+1);
% Set up matrix for recursions, initialise as ARG terminates at nl=1, 
%and need R recombs to have occurred by this point, none of which are open.
p_matrix(1, 1, R+1) = 1;

%The galled prob reutrn function returns the probability if it exists (the
%indices make physical sense) and 0 if else.
for r = fliplr(0: R)
    for r0 = 0:r
        for nl = 2:n+r0
            p_matrix(nl, r0+1, r+1) =(r0*galled_prob_return(p_matrix, nl-1, r0-1, r, n, R) +...
                ((nl-2*r0)*(nl-1-2*r0)/2 + (nl - 2*r0)*2*r0)* galled_prob_return(p_matrix, nl-1, r0, r, n, R) + ...
                rho/2*(nl - 2*r0)*galled_prob_return(p_matrix, nl+1, r0+1, r+1, n, R)) / (nl/2*(nl-1+rho));
        end
    end
end

% Final probability given by p_matrix(n, r0=0, r=0)
p_sol = p_matrix(n, 1, 1);
p_matrix_sol = p_matrix;
end


