function [q_sol, q_matrix_sol] = solve_number_recombs(n, R, rho)
% Solves the recursion for the probabilities that a total of R
% recombinations occur in a sample history size n.
% Returns single probability and also whole matrix calculated in
% recursions.
% q(nl, r) is the probability conditioned on r recombs happened so far
% in history, with nl lineages remaining

q_matrix = zeros(n+R, R+1);
% Set up matrix for recursions, initialise as ARG terminates at nl=1, 
%and need R recombs to have occurred by this point.
q_matrix(1, R+1) = 1;

%Solve for when all recombs have occurred.
for nl = 2:n+R
    q_matrix(nl, R+1) = (nl - 1)/(nl - 1 + rho)*q_matrix(nl - 1, R+1);
end

% Solve for rest of matrix.
for r = fliplr(0: R-1)
    for nl = 2:n+r
        q_matrix(nl, r+1) = 1/(nl- 1 + rho) * ((nl-1)*q_matrix(nl-1, r+1) + rho*q_matrix(nl+1, r+2));
    end
end

% Final probability given by q_matrix(n, 0)
q_sol = q_matrix(n, 1);
q_matrix_sol = q_matrix;
end


