
load('DataDuroturning_Ly10h1over16_tf32_eee3_resist0.mat');
vv = [0,0];  h = 1/16;
% Comment: In this example, I use i = 8,66,115,159,196 to generate the Duroturning behavior for old BC in our paper.  
i = 8;
figure(i);
t = 0.1*i; pretitle  = sprintf('time t = %.4f', t);     
% load the snapshot data
Phi  = Phisnap(:,:,i);   Xbar = xbars(i,1); Ybar = ybars(i,1);
tmp1 = dragXsnap(:,:,i); tmp2 = dragYsnap(:,:,i);
% plot the result
contour(X,Y,Phi0,vv,'w');
hold on;
contour(X,Y,Phi,vv,'g');
scatter(Xbar,Ybar);
leg = 1/(2*h); 
quiver(X(1:leg:end,1:leg:end),Y(1:leg:end,1:leg:end),tmp1(1:leg:end,1:leg:end)/2,tmp2(1:leg:end,1:leg:end)/2); 

axis([-11 11 -5 10]); 
pbaspect([22 15 1]);
set(gca,'XTick',[],'XColor','w','YTick',[],'YColor','w');
set(gca,'position',[0,0,1,1]);