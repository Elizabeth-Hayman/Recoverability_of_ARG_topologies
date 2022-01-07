function P = solve_full_topology(R,n,N0,rho,theta,z)
%Solve for the ful1 ARG topologies.

%Manipulate theta into helpful forms
t=theta/2;
ta=t*z;
tb=t*(1-z);

REF=ref_create(R);
length_ref = length(REF(1,:));
R1=R+1;
p=zeros(N0+1, n+R, n+R1, R1, R1, length_ref,length_ref); %p(no,nl,nf,r0,r, lefthand, righthand)

% LEFT i - recomb occurs,  j- first coalesc, k1- associated mutation, k2 - non-associated mutation, l1-
% second coalesce, m1- non-associated mutation, m2- associated mutation,
% l2- third mutation.
% e- open lin, a - first coalesc, b- associated mutation, bb non-associated mutation c-
% second coalesce, d- non-associated mutation, dd associated mutation

p(:, 1, 1, 1, R+1, 1, 1) = ones(1, N0+1); %nl=1, nf=0, r0=0, r=R, states=0
p(:, 1, 2, 1, R+1, 1, 1) = ones(1, N0+1); %nl=1, nf=1, r0=0, r=R, states=0 for completeness,
%though termination at n=1 means no further mutations would arise. 

for n0=0:N0
    for r = fliplr(0: R)
    for r0 = 0:r
        for nl=2:n+r0
        for nf=fliplr(0:nl-2*r0)
            %Restrictions on range on nf, nl due to number of recombination
            %loops open. Note nf does NOT include any recombination
            %lineages.
            nu = (nl-2*r0-nf); % number of non-recombinant lineages not ready to coalesce
            for i=0:r0
            for j=0:r0
            for k2=0:r0
            for k1=0:r0
            for l2=0:r0
            for l1=0:r0
            for m2=0:r0
            for m1=0:r0
            for e=0:r0
            for a=0:r0
            for b2=0:r0
            for b1=0:r0
            for c1=0:r0
            for c2=0:r0
            for d2=0:r0
            for d1=0:r0
                if i+j+k1+k2+l1+l2+m1+m2 == r0 && e+a+b1+b2+c1+c2+d1+d2 == r0
                    %only continue if we have a valid tuple of recombinant
                    %states
                    ref1=ref_return(REF,i,j,k1,k2,l1,m1,m2,l2); %ref mat finds the index for the state of the ARG
                    ref2=ref_return(REF,e,a,b1,b2,c1,d1,d2,c2); %indentical and independent to find the index for right hand lineages
                    
                    Rate = (nl*(nl-1)/2+rho/2*nl+t*(nl-nf-2*r0+c2+l2)+ta*(j+k2+l1+a+c1+d2)+tb*(j+l1+m2+a+b2+c1));
                    
                    %coalesce of non-recomb lins
                    Coal_NR = nf*(nf-1)/2*top_prob_return(p,n0,nl-1,nf-2,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + ...
                              nf*nu*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + ...
                              nu*(nu-1)/2*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF);
                    
                    %new recombination
                    Recomb = rho/2*(nu*top_prob_return(p,n0,nl+1,nf,r0+1,r+1,i+1,j,k1,k2,l1,m1,m2,l2,e+1,a,b1,b2,c1,d1,d2,c2,n,R,REF)+... %new gal1ed recomb
                                    nf*top_prob_return(p,n0,nl+1,nf-1,r0+1,r+1,i+1,j,k1,k2,l1,m1,m2,l2,e+1,a,b1,b2,c1,d1,d2,c2,n,R,REF)); %need to separate cases
                    
                    %mutation of a non-recomb lineage    
                    Mut_NR = t*nu*top_prob_return(p,n0,nl,nf+1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF);
                    
                    %coalesce of one recomb and one non-recomb lin
                    Coal_R = i*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i-1,j+1,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i-1,j+1,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ...
                        j*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ...
                        e*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e-1,a+1,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e-1,a+1,b1,b2,c1,d1,d2,c2,n,R,REF))+...
                        a*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ... 
                        k1*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1-1,k2,l1+1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1-1,k2,l1+1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) +... 
                        b1*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1-1,b2,c1+1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1-1,b2,c1+1,d1,d2,c2,n,R,REF)) +...
                        k2*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j+1,k1,k2-1,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j+1,k1,k2-1,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) +...
                        b2*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a+1,b1,b2-1,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a+1,b1,b2-1,c1,d1,d2,c2,n,R,REF)) +...
                        m1*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1-1,m2,l2+1,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1-1,m2,l2+1,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) +... 
                        c1*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ...
                        m2*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1+1,m1,m2-1,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1+1,m1,m2-1,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) +...
                        c2*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ...
                        l1*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF))+...
                        d1*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2+1,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2+1,n,R,REF)) + ...
                        l2*(nf*top_prob_return(p,n0-1,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + nu*top_prob_return(p,n0-2,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF)) + ...
                        d2*(nf*top_prob_return(p,n0,nl-1,nf-1,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1+1,d1,d2-1,c2,n,R,REF) + nu*top_prob_return(p,n0-1,nl-1,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1+1,d1,d2-1,c2,n,R,REF));
                    
                    %mutation of recomb lins    
                    Mut_R = ta*j*top_prob_return(p,n0,nl,nf,r0,r,i,j-1,k1+1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + tb*a*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a-1,b1+1,b2,c1,d1,d2,c2,n,R,REF) + ... 
                        ta*k2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1+1,k2-1,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + tb*b2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1+1,b2-1,c1,d1,d2,c2,n,R,REF) + ...
                        tb*j*top_prob_return(p,n0,nl,nf,r0,r,i,j-1,k1,k2+1,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + ta*a*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a-1,b1,b2+1,c1,d1,d2,c2,n,R,REF) + ...
                        tb*l1*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1-1,m1+1,m2,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + ta*c1*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1-1,d1+1,d2,c2,n,R,REF) + ...
                        tb*m2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1+1,m2-1,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + ta*d2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1+1,d2-1,c2,n,R,REF) + ...
                        ta*l1*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1-1,m1,m2+1,l2,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + tb*c1*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1-1,d1,d2+1,c2,n,R,REF) + ...
                        t*l2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1+1,m2,l2-1,e,a,b1,b2,c1,d1,d2,c2,n,R,REF) + t*c2*top_prob_return(p,n0,nl,nf,r0,r,i,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1+1,d2,c2-1,n,R,REF);
                    
                    
                    %recoalescence of a recombination loop. If statement
                    %prevent errors being thrown when r0=0
                    if r0 == 0
                        Coal_RR = 0;
                    else
                        Coal_RR = m1/r0*(e*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e-1,a,b1,b2,c1,d1,d2,c2,n,R,REF) + b1*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1-1,b2,c1,d1,d2,c2,n,R,REF) +...
                                b2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1,b2-1,c1,d1,d2,c2,n,R,REF) +  d2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1,b2,c1,d1,d2-1,c2,n,R,REF) +...
                                 d1*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + c1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1,b2,c1-1,d1,d2,c2,n,R,REF) +...
                                 c2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + a*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1-1,m2,l2,e,a-1,b1,b2,c1,d1,d2,c2,n,R,REF)) +...
                            l2/r0*(e*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e-1,a,b1,b2,c1,d1,d2,c2,n,R,REF) + b1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1-1,b2,c1,d1,d2,c2,n,R,REF) +...
                                b2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2-1,c1,d1,d2,c2,n,R,REF) +  d2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2,c1,d1,d2-1,c2,n,R,REF) +...
                                 d1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + c1*top_prob_return(p,n0-2,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2,c1-1,d1,d2,c2,n,R,REF) +...
                                 c2*top_prob_return(p,n0-2,nl-1,nf-2,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + a*top_prob_return(p,n0-2,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a-1,b1,b2,c1,d1,d2,c2,n,R,REF)) +...
                            d1/r0*(i*top_prob_return(p,n0,nl-1,nf,r0-1,r,i-1,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + k1*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1-1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + ...
                                k2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2-1,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + m2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2-1,l1,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) +...
                                j*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j-1,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF) + l1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2,l2-1,e,a,b1,b2,c1,d1-1,d2,c2,n,R,REF)) +...
                            c2/r0*(i*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i-1,j,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + k1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1-1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + ...
                                k2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2-1,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + m2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2-1,l1,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) +...
                                j*top_prob_return(p,n0-1,nl-1,nf-1,r0-1,r,i,j-1,k1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF) + l1*top_prob_return(p,n0,nl-1,nf-2,r0-1,r,i,j,k1,k2,l1-1,m1,m2,l2,e,a,b1,b2,c1,d1,d2,c2-1,n,R,REF)) +...
                            k1/r0*(b1*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1-1,k2,l1,m1,m2,l2,e,a,b1-1,b2,c1,d1,d2,c2,n,R,REF) + d2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1-1,k2,l1,m1,m2,l2,e,a,b1,b2,c1,d1,d2-1,c2,n,R,REF) +...
                                c1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1-1,k2,l1,m1,m2,l2,e,a,b1,b2,c1-1,d1,d2,c2,n,R,REF))+...
                            m2/r0*(b1*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2-1,l1,e,a,b1-1,b2,c1,d1,d2,c2,n,R,REF) + d2*top_prob_return(p,n0,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2-1,l1,e,a,b1,b2,c1,d1,d2-1,c2,n,R,REF) +...
                                c1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1,m1,m2-1,l1,e,a,b1,b2,c1-1,d1,d2,c2,n,R,REF))+...
                            l1/r0*(b1*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1-1,m1,m2,l2,e,a,b1-1,b2,c1,d1,d2,c2,n,R,REF) + d2*top_prob_return(p,n0-1,nl-1,nf,r0-1,r,i,j,k1,k2,l1-1,m1,m2,l2,e,a,b1,b2,c1,d1,d2-1,c2,n,R,REF) +...
                                c1*top_prob_return(p,n0-2,nl-1,nf,r0-1,r,i,j,k1,k2,l1-1,m1,m2,l2,e,a,b1,b2,c1-1,d1,d2,c2,n,R,REF));
                    end
                        
                        
                        p(n0+1, nl, nf+1, r0+1, r+1, ref1, ref2) = (Coal_NR + Coal_R + Coal_RR + Mut_R + Mut_NR + Recomb)/Rate;
                end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end        
        end                                         
        end
    end
    end
end
P=p;
%P(N0+1, n, n+1, 0, 0, 1, 1); gives the prob starting from a sample size n
end