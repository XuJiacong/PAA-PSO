function [ ind ] = fuzzyclassification( f,ind_1) 
queen=[1 2 3 4]; 
state=[s1(f),s2(f),s3(f),s4(f)]; 
if (length(find(state~=0))==1) 
    [~, ind] = max(state); 
else 
    record=find(state>0); 
    if record(1)==ind_1 
        ind=ind_1; 
    elseif record(2)==ind_1 
        ind=ind_1; 
    else 
        for i=1:1:4 
            j=mod(ind_1+i,4); 
            if j==0 
                j=4; 
            end 
            if queen(j)==record(1) 
                ind=record(1); 
                break
            elseif queen(j)==record(2) 
                 ind=record(2); 
                 break
            end 
        end 
    end 
end 
end 
 