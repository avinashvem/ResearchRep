function output=repeat_values(arr)
arr=reshape(arr,1,[]);
output=[];
if length(unique(arr))~=length(arr)
    sort_arr=sort(arr);
    diff=[sort_arr sort_arr(end)+1]-[0 sort_arr];
    zero_idx=find(diff==0);
    for i=1:numel(zero_idx)
        idx=zero_idx(i);
        frq=2;
       while diff(idx+1)==0
           idx=idx+1;
           frq=frq+1;
       end
       output=[output;[sort_arr(idx) frq]]; 
    end
end