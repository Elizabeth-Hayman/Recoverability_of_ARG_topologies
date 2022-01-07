function prob=top_prob_return(prob_matrx,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,ref) %short function outputs 
%Returns non-trivial probabilities only if the indices make physical sense

nl_min = max(1, 2*r0);
if nl_min<=nl && nl<=n+r0 && 0<=nf && nf<=nl && 0<=n0 && 0<=r0 && r0<=r && r<=R && ...
        i+j+k1+k2+l1+l2+m1+m2 == r0 && e+a+b1+b2+c1+c2+d1+d2 == r0 && ...
        0<=i && 0<=j && 0<=k1 && 0<=k2 &&  0<=l1 && 0<=m1 && 0<=m2 && 0<=l2 && ...
        0<=e && 0<=a && 0<=b1 && 0<=b2 && 0<=c1 && 0<=d1 && 0<=d2 && 0<=c2

    ref1=ref_return(ref,i,j,k1,k2,l1,m1,m2,l2);
    ref2=ref_return(ref,e,a,b1,b2,c1,d1,d2,c2);
    
    if isempty(ref1)
        [i,j,k1,k2,l1,m1,m2,l2]
    end
    if isempty(ref2)
        [e,a,b1,b2,c1,d1,d2,c2]
    end
    
    p=prob_matrx(n0+1,nl,nf+1,r0+1,r+1,ref1,ref2);
                                   
else
    p=0;
end

prob=p;
end