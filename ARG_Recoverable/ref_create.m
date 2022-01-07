function REF=ref_create(R)
%Creates a 2D array with assigning each set of states the left recombinant
%lineages are in to an integer. Restricted by up to R total recombinations
%allowed to be open any one time.

%Code works identically with right hand recombinant lineages as these
%behave identically and independently to the left hand.

%Output will look like:
% 0   1   10   
% 1   2    3

%where the first two entries should be enterpreted as 00  01  10

size=[9,45,165,495,1287,3003,6435];
if R>0
    length_ref=size(R);
else
    length_ref=1;
end
ref=zeros(2,length_ref);
ref(2,:)=linspace(1,length_ref,length_ref);
count=1;

    for i=0:R
        for j=(0:R)
            for kk=(0:R)
                for k=(0:R)
                    for l=0:R
                        for mm=(0:R)
                            for m=(0:R)
                                for ll=0:R
                                if i+j+k+kk+m+mm+l+ll<=R
                                    ref(1,count)=string(i)+string(j)+string(k)+string(kk)+string(l)+string(m)+string(mm)+string(ll);
                                    count=count+1;
                                end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    REF=ref;
end