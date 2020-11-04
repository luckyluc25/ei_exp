function a=select_action(p)

rng('shuffle')

if unique(p)<length(p)
    flag=1;
end

p_sum=cumsum(p);
rand_num=rand(1);

interval=find(rand_num>p_sum);
if isempty(interval)
    a=1;
else
    a=interval(length(interval))+1;
end

end
