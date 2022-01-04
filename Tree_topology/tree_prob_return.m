function p = tree_prob_return(p_matrix, n0, nl, nf, n, N0)
% Function to return the value of an array if position called exists, and
% to return 0 if else.

if 0<=nf && nf<=nl && nl<=n && 0<=n0 && n0<=N0
    p_value = p_matrix(n0+1, nl, nf+1);
    
else
    p_value = 0;
end
p=p_value;
end