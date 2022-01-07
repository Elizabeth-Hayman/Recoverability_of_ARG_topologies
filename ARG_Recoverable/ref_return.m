function index=ref_return(ref_matrix,I,J,K1,K2,L1,M1,M2,L2)
% Takes a set of left hand states, turns it into a reference string and
% retrieves the index from a given reference matrix.

reference=str2double(string(I)+string(J)+string(K1)+string(K2)+string(L1)+string(M1)+string(M2)+string(L2));
    
ind=find(ref_matrix(1,:)==reference,1);
    
index=ref_matrix(2,ind);
end 