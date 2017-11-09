cb = repmat([0 0 0],length(D),1);
scatter(D(1,:),D(2,:),40,cb);title('Small supervise');
hold on
cs2 = repmat([0 1 0],5,1);
cs1 = repmat([1 0 0],5,1);
scatter(D(1,n(1:4:20)),D(2,n(1:5)),80,cs1,'filled','d');
scatter(D(1,n(461:4:480)),D(2,n(1:5)),80,cs2,'filled','d');
figure
c = FD == 1;
c1 = repmat([0 1 0],247,1);
c3 = repmat([1 0 0],length(D),1);
c3(c,:) = c1;
scatter(D(1,:),D(2,:),40,c3);title('Ground Truth');
