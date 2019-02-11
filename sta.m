function shape = sta(rho,stim,dt,tau)
%   rho--sequence of events/nonevents specified by 0/1
%   stim--stimulus
%   tau--an array of time window sizes
tau=tau/dt; %convert window size from time to points
shape=cell(length(tau),1);
for i=1:length(tau)
    shape{i}=zeros(tau(i),size(stim,2)); %initialize shape for each window size
end
count_mat=zeros(length(tau),1); %record the number of spikes events for each shape{i}
for i=1:length(rho)
    if rho(i)==1
        for j=1:length(tau)
            if i-tau(j)>=1 %only pick events with enough points preceding 
                shape{j}=shape{j}+stim(i-tau(j):i-1,:); % for each time window size, add stimulus to the corresponding cell in shape
                count_mat(j)=count_mat(j)+1;
            end
        end
    end
end

for i=1:length(tau)
    shape{i}=shape{i}/count_mat(i); %average over the window size
end

