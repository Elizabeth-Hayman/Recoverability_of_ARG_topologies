function p = galled_prob_return(p_matrix, nl, r0, r, n, R)
% Function to return the value of an array if position called exists, and
% to return 0 if else.

if 2*r0<=nl && nl<=n+r0 && 0<=r0 && r0<=r && r<=R
    p_value = p_matrix(nl, r0+1, r+1);
    
else
    p_value = 0;
end
p=p_value;
end