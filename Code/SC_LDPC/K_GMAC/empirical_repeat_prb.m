function output_prb=empirical_repeat_prb(Ka,M,numSims)

Mat_coll=zeros(1,5);
for i=1:numSims
X_coeffs=sort(randi(M,[1,Ka]));
repeatV=repeat_values(X_coeffs);
uniq_coeffs=unique(X_coeffs);

switch length(uniq_coeffs)
            case Ka
                collision_flag=0;
            case Ka-1
                collision_flag=1;
            case Ka-2
                if length(repeatV)==1
                    collision_flag=3;
                elseif length(repeatV)==2
                    collision_flag=2;
                else
                    collision_flag=-1;
                end
    otherwise
        collision_flag=-1;
end
                    
    if collision_flag==-1
         Mat_coll(5)=Mat_coll(4)+1;
    else
         Mat_coll(collision_flag+1)=Mat_coll(collision_flag+1)+1;
    end
end
output_prb=Mat_coll./numSims;
end