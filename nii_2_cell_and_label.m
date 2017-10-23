filepath='F:\ADNI\ALL';
filelist=dir(filepath);%All nii data without scaled_2,there are 818 files
N=length(filelist);

data=cell(N,1);
for i=1:N
    filename=filelist(i);
    if length(filename.name)>2  %maybe there is system file name '.'
        V=spm_vol(fullfile(filepath,filename.name));
        data{i,1}=spm_read_vols(V);
    end
end
data(cellfun(@isempty,data))=[];

%{
V1=spm_vol(path1);x1=spm_read_vols(V1);
V2=spm_vol(path2);x2=spm_read_vols(V2);
V3=spm_vol(path3);x3=spm_read_vols(V3);
V4=spm_vol(path4);x4=spm_read_vols(V4);
y1=x1(:);y1=y1/max(y1);%NL
y2=x2(:);y2=y2/max(y2);%AD
y3=x3(:);y3=y3/max(y3);%MCI
y4=x4(:);y4=y4/max(y4);%NL
order=1;
level=1;
[H]=spfilter(order);
[fcoef1]=framedec3d(x1,H,level);
[fcoef2]=framedec3d(x2,H,level);
[fcoef3]=framedec3d(x3,H,level);
[fcoef4]=framedec3d(x4,H,level);
data=zeros(27,3);
for i=1:27
    data(i,1)=sum(sum(sum(fcoef1{i}.^2)));
    data(i,2)=sum(sum(sum(fcoef2{i}.^2)));
    data(i,3)=sum(sum(sum(fcoef3{i}.^2)));
    %data(i,4)=sum(sum(sum(fcoef4{i}.^2)));
end
%}