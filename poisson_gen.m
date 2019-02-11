function [spikes, intervals] = poisson_gen(Ttot,dt,r_est,tau)
count=0;
r=r_est;
intervals=[];
spikes=[];

if length(r_est)==1
    for i=1:Ttot/dt
        if tau~=0  %update r at each time step
            r=r+(r_est-r)/tau*dt;
        end
        if r*dt>rand
            count=count+1;
            spikes(count)=i*dt;
            if count>1
                intervals(count-1)=spikes(count)-spikes(count-1);
            end
            if tau~=0 %update r after spike
                r=0;
            end
        end
    
    end
else  %if r_est is an array, i.e. a function of time
    num_r=length(r_est); %get the number of different r, 
    % each will be used in an interval further divided by dt in the loop
    dT=Ttot/num_r; % the number of bins in calculating r as a function of time
    for i=1:num_r
        for j=1:dT/dt
            
            if r_est(i)*dt>rand % unit /data point for this question
                count=count+1;
                spikes(count)=(i-1)*dT+j*dt; %the spike time is the jth dt of the ith bin
            end
        end
    end
    % we don't care about the intervals in this case
end


end