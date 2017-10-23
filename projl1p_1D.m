%orthogonal projection of one or more vectors onto nonnegative face of l1 a-ball
%sz = matrix with columns to be projected
%a = ball radius
%len = vector length
%vecs = number of vectors
%version 3: removed superfluous sort
% u = projl1p_1D(u,1); call in SplitBregGraphClassK
function z = projl1p_1D(sz,a) 
[ny,len] = size(sz);
vecs=ny;
sz=sz';
%length 1 special case
if (len==1)
   z = a*ones(ny,1);
   return
end
nx=1;
%initializations
C = zeros(1,vecs);
rind = (1:len)'; %column vector of row indices
cind = (1:ny*nx); %row vector of column indices
Row = rind(rind,ones(1,vecs)); %matrix of row indices
Col = cind(ones(len,1),cind); %matrix of column indices
thr = zeros(1,vecs); %row vector of column thresholds
z = zeros(len,vecs);

%sort down the columns
[sz,SI] = sort(sz,1); %also keep track of sorting indices
SI = SI + len*(Col-1);

%build cumsum matrix (value of sum from current row to end in that column)
S = cumsum(sz(end:-1:1,:));
S = S(end:-1:1,:);

%set km and kM (these are actually row indices, different for each column)
km = ones(1,vecs);
kM = len*ones(1,vecs);
k = km;

%either km is strict lower bound or column is done
%project onto planes at km=1 threshold (all columns active)
thr = (S(1,:)-a)/len;
p = sz - thr(ones(len,1),cind);

%if >= 0 then column no longer active
Y = min(p,[],1)<0;

%if column done already then set kM = km
kM(~Y) = km(~Y);

%bisect to find threshold k
act = sum(Y);
while ( act > 0 )
   %bisect (ceil to integer)
   k(Y) = ceil((kM(Y)+km(Y))/2);
   
   %get indices
   ind = len*(cind(Y)-1)+k(Y);
   
   %check condition for each active column
   C(Y) = sz(ind)-(S(ind)-a)./(len-k(Y)+1)<0;
   
   %if negatives, increase km, othersise decrease kM       
   km(Y&C) = k(Y&C);
   kM(Y&~C) = k(Y&~C);
   
   %column done if kM = km + 1 (update Y and act)
   Y(Y) = kM(Y)-km(Y)-1~=0;
   act = sum(Y);
end

%apply threshold kM to sz
ind = len*(cind-1)+kM;
thr = (S(ind)-a)./(len-kM+1);
p = sz - thr(ones(len,1),cind);
p = p.*(Row >= kM(ones(len,1),cind)); %threshold

%return z
z(SI) = p;
z=reshape(z',ny,len);