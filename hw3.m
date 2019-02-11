%% Q2
clear all
load('fly_data.mat')
rl=length(rho)/2;
sl=length(stim)/2;
window=[10/500,30/500,50/500,100/500,200/500];
dt=1/500
q2shape=sta(rho(1:rl),stim(1:sl),dt,window);
figure
hold on
xlabel('time (s)')
ylabel('stimulus')
wl=length(window)
for i=1:wl
    plot(linspace(window(wl)-window(wl-i+1),window(wl),window(wl-i+1)/dt),q2shape{wl-i+1},'DisplayName','WindowSize');
    
end
hold off

%from the graph it is clear that 100/500 would be a better time window,
%because it contains the entire process of change of the stimulus, whereas
%a smaller or a bigger one would either miss a part or contains unnecessary
%parts.

%% Q3
sta_shape=q2shape{4}; %since 100/500s is a better time window
%calculating P(L)
L=zeros(length(stim)/2-length(sta_shape)+1,1); %initialize the overlap array
for i=1:length(L)-length(sta_shape)+1 %while the number of data points left in stim is >= the size of the filter
    sum=0;
    for j=1:length(sta_shape)
        sum=sum+sta_shape(j)*stim(i-1+j); %dot product between a moving filter and the stimulus
    end
    L(i)=sum;
end

L=L/var(stim(1:length(stim)/2)); %normalize the overlap
[P_L,edge_L]=histcounts(L,'BinWidth',1,'Normalization','probability');
figure
histogram(L,edge_L,'Normalization','probability')
title('P(L)')
%calculating P(L|spike)
LgivenSp=[];
count=1;
for i=1:length(rho)/2
    if rho(i)==1
        if i>length(sta_shape)  %ignore cases where there's no enough data points before a spike in stim
            
            LgivenSp(count)=sta_shape'*stim(i-length(sta_shape):i-1); %vector dot product
            count=count+1;
        end
    end
end

LgivenSp=LgivenSp/var(stim(1:length(stim)/2));


[P_LgivenSp,edge_LgivenSp]=histcounts(LgivenSp,edge_L,'Normalization','probability'); %use the same edges as L
figure
histogram(LgivenSp,edge_L,'Normalization','probability');
title('P(L|Spike)')
P_SPgivenL=P_LgivenSp./P_L;
%calculating P(spike)
f=find(rho(1:length(rho)/2)==1);
firing_rate=length(f)/(length(rho)/2);
P_Sp=firing_rate*exp(-firing_rate);

%P(spike|L)
P_SPgivenL=P_Sp*P_LgivenSp./P_L;
figure
bar(edge_L(length(edge_L)-length(P_SPgivenL)+1:length(edge_L)),P_SPgivenL); %adjust edge_L to make dimension match; effectively this abandon the leftmost edge
title('P(spike|L)')

%% Q4
%generate instantaneous firing rate for window_size=100 data points, using P(spike|L)
v=var(stim(1:length(stim)/2));
test_stim=stim(length(stim)/2+1:length(stim)); %extract the second half of stim
r_est=[];
for i=1:length(test_stim)-length(sta_shape)+1 %again, not counting the end where the number of data points left is less than the size of the filter
    if i>length(sta_shape)
        l=sta_shape'*test_stim(i-length(sta_shape):i-1)/v; %calculate the overlap
        
        r_est(i)=P_SPgivenL(ceil(l)-edge_L(1)); %ceil find out which bin l is in; the convert it to the corresponding index in P_SPgivenL
    end
end
[spikes,intvl]=poisson_gen(length(r_est),1,r_est,0);

%visualize; rasterplot for every 1000 points
test_rho=rho(length(rho)/2+1:length(rho));

min=1;
max=1; % for each 1000 points in time, max is the index of the last element in spikes that is in the time internal, and min the first
% generate 10 plots


for i=1:10
    if max>length(spikes) % making sure index does not exceed dimension
       break
    end
    
    while (i-1)*1000<spikes(max+1)&&spikes(max+1)<=i*1000
        max=max+1;
        if spikes(max)>i*1000
            max=max-1;  %the last max already exceeds the interval and therefore needs to be subtracted
            break
        end
    end
   
    model=2*ones(max-min+1,1);
    figure
    hold on 
    scatter(spikes(min:max),model);
    min=max+1; %the one after the max for the previous interval is the min for the next
    scatter([(i-1)*1000+1:i*1000],test_rho((i-1)*1000+1:i*1000));  
    title(['spike train: model vs data ',num2str(i)])
    xlabel('time points (2ms)')
    ylim([0.1,2])
    legend('model','data')
    hold off
    
end

%% Q5
clear all
load('533l006.p05_stc.mat')
stim5=stim;
spikes_per_frm5=spikes_per_frm;
load('543l021.p07_stc.mat')
stim7=stim;
spikes_per_frm7=spikes_per_frm;
%spike triggered average
temp5=sta(spikes_per_frm5,stim5,1,30); %time window of 100 points
sta5=temp5{1}; %the unnecesary complication arises because my sta produces a cell array
temp7=sta(spikes_per_frm7,stim7,1,30);
sta7=temp7{1};
figure
imagesc(sta5)
colormap gray
colorbar
xlabel('spatial dimensions')
ylabel('time')
title('STA 05')
figure
imagesc(sta7)
colormap gray
colorbar
xlabel('spatial dimensions')
ylabel('time')
title('STA 07')
%Spike triggered covariance
[stc5,V5,D5]=stc(stim5,spikes_per_frm5,sta5);
[stc7,V7,D7]=stc(stim7,spikes_per_frm7,sta7);
%plot eigenvalue spectrum
figure
hold off
scatter([1:length(D5)],D5);
title('Eigenvalue Spectrum 05');
xlabel('Eigenvalue number');
ylabel('variance');
figure
scatter([1:length(D7)],D7);
title('Eigenvalue Spectrum 07');
xlabel('Eigenvalue number');
ylabel('variance');
%There is one significant eigenvalue for P05 and six for P07.  
for i=[1,2,length(V5)] %putting 2 here just to show that the second is not significant
    v5=reshape(V5(:,i),[size(sta5,1),size(sta5,2)]); %convert vector into matrix
    figure
    imagesc(v5)
    colormap gray
    colorbar
    xlabel('spatial dimensions')
    ylabel('time')
    title(['Eigenvector',num2str(i),' 05'])
end
for i=[1:6,length(V7)]
    v7=reshape(V7(:,i),[size(sta7,1),size(sta7,2)]);
    figure
    imagesc(v7)
    colormap gray
    colorbar
    xlabel('spatial dimensions')
    ylabel('time')
    title(['Eigenvector',num2str(i),' 07'])
end    
%For P05 since there is only one significant eigenvalue, the spike
%triggered covariance does not offer much more infromation than STA.
%However for P07 the spike triggered average shows only a horizontal dark
%bar, while the STC shows how it is composed of different filters with
%richer internal structure (light/dark bar in between two dark/light shades).


