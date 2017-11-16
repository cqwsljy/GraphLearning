% »­mnistÊı¾İÍ¼

load mnistAll.mat

img0 = zeros(28*10,28*10);
for k = 0:9
    index = find(FD == k);
    index =  index(randperm(length(index)));
    for j = 1:10
       data = D(:,index(j));
       data = reshape(data,28,28)';
       img0(28*(j-1)+1:28*j,28*(k)+1:28*(k+1)) = data;
    end
end