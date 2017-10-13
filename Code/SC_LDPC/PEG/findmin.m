function [num,dc_min_id]=findmin(vec,deg)
m=length(vec);
dc_min_id=zeros(1,m);
num=0;
for i=1:m
    if vec(i)~=0
        break
    end
end
if vec(i)>0    
    dc_min_id=vec(i);
    dc_min=deg(vec(i));


    for t=i+1:m
        if vec(t)>0
            if deg(vec(t))<dc_min
                dc_min=deg(vec(t));
%                 dc_min_id=vec(t);
            end
        end
    end
end
for i=1:m
    if vec(i)>0
        if deg(vec(i))==dc_min
            num=num+1;
            dc_min_id(num)=vec(i);
        end
    end
end
        
        
end




