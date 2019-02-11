function [cov_mat,V,D]= stc(stim,spikes,sta)
%return the spike triggered covariance matrix
%stim--stimulus
%spikes--binary spike train
%sta--spike triggered average

dim=size(sta,1)*size(sta,2);
cov_mat=zeros(dim);
count=0; %keeping track of the total number of summation N
sta_vec=reshape(sta,[dim,1]); %convert the matrix into a column vector
for i=1:length(spikes)
    if spikes(i)==1
        if i>length(sta)
            stim_vec=reshape(stim(i-length(sta):i-1,:),[dim,1]); %convert the matrix into a column vector
            cov_mat=cov_mat+(stim_vec-sta_vec)*(stim_vec-sta_vec)';
            count=count+1;
        end
    end
end
cov_mat=1/(count-1)*cov_mat;            

[V,D]=eig(cov_mat);  %the sorting code is from matlab documentation
[D,ind]=sort(diag(D),'descend'); %https://www.mathworks.com/help/matlab/ref/eig.html

V=V(:,ind);





end

