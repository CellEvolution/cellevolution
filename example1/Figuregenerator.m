load('DataFragmentpush_tf20_tss200_h1over16.mat')
vv = [0,0];  h = 1/16;
% Comment: In this example, I use i = 1,4,16,26, 70 to generate the Push behavior. 

i =26; figure(i);
t = 0.1*i; pretitle  = sprintf('time t = %.4f', t);     
% load the snapshot data
Phi  = Phisnap(:,:,i);   Xbar = xbars(i,1); Ybar = ybars(i,1);
tmp1 = dragXsnap(:,:,i); tmp2 = dragYsnap(:,:,i);
% plot the result
contour(-X,Y,Phi0,vv,'r');
hold on;
contour(-X,Y,Phi,vv,'g');
scatter(-Xbar,Ybar);
leg = 1/(2*h); 
quiver(-X(1:leg:end,1:leg:end),Y(1:leg:end,1:leg:end),-tmp1(1:leg:end,1:leg:end)/2,tmp2(1:leg:end,1:leg:end)/2); 
axis([-4 8 -4 4]); 
pbaspect([12 8 1]);

figure(i)
set(gca,'XTick',[],'XColor','w','YTick',[],'YColor','w');
set(gca,'position',[0,0,1,1])