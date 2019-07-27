% This code extends the velocity define on the interface to whole computational domain
% Use the formula  v=n*Sn-<n*Sn>_{\partial \Omega(t)} in whole domain
% zhiwen zhang/20141121
% Caution: obtain the length of the interface and average via the regularized delta function!
function [Vn,NormalVecX,NormalVecY,nSn,nSn_ave,Ux,Uy,Vx,Vy,tp1,tp2] = ...
         RadialModelTwo_NormalStressV6CORRECTEDoffset(X,Y,Phi0,Phi,U,V,p,mu,lmd,theta,beta,zz,zz1,eee,offset,lambda,K,gamma,resist,expN,Hev,Del,Xbart,Ybart)
%% parameter
alpha     = 2*theta*(mu+lmd)*0;
%beta      = 1;
h         = (p(2,2)-p(1,2));  
Xbar      = sum(sum(Hev(Phi).*X))/sum(sum(Hev(Phi)));   % centroid of the cell
Ybar      = sum(sum(Hev(Phi).*Y))/sum(sum(Hev(Phi)));
% data structure to record the value and distance for each nodal point (x_i,y_j), Right,Left,Up and Bottom
UR = U(:,[2:end 1]);       VR = V(:,[2:end 1]);            % For x-direction
UL = U(:,[end 1:end-1]);   VL = V(:,[end 1:end-1]);        % For x-direction
UU = U([2:end 1],:);       VU = V([2:end 1],:);            % For y-direction
UB = U([end 1:end-1],:);   VB = V([end 1:end-1],:);        % For y-direction
%% compute
%%%%%%%%%% calculate the 1st order partial derivatives Ux,Uy,Vx and Vy %%%%%%%%%%
Ux = (UR - UL)/2/h;      Uy = (UU - UB)/2/h;  
Vx = (VR - VL)/2/h;      Vy = (VU - VB)/2/h; 
%%%%%%%%%% special treat four boundary sides derivatives 
Ux(:,1)=(U(:,2)-U(:,1))/h; Ux(:,end) = (U(:,end)-U(:,end-1))/h;     
Vx(:,1)=(V(:,2)-V(:,1))/h; Vx(:,end) = (V(:,end)-V(:,end-1))/h;     
Uy(1,:)=(U(2,:)-U(1,:))/h; Uy(end,:) = (U(end,:)-U(end-1,:))/h;      
Vy(1,:)=(V(2,:)-V(1,:))/h; Vy(end,:) = (V(end,:)-V(end-1,:))/h;    
%%%%%%%%% Heviside and Delta of the level set function %%%%%%%%%%%
H_Phi  = Hev(Phi); 
D_Phi  = Del(Phi);
%%%%%%%%% calculate the stress tensor %%%%%%%%
sigma11 = (lmd+2*mu)*Ux + lmd*Vy       + alpha*H_Phi;
sigma22 =  lmd*Ux       +(lmd+2*mu)*Vy + alpha*H_Phi;
sigma12 =  mu*(Uy+Vx);
%% compute
%%%%%%%%%% calculate the velocity on the interface and label the displacement value on the irregular grid %%%%%%%%%% 
% calculate the radial direction on each nodal point %%%%%%%%%%
erx0 = X - Xbar;        ery0 = Y - Ybar; 
RadialVec = sqrt(erx0.^2 + ery0.^2 + 1e-16); 
erx = erx0./RadialVec;  ery = ery0./RadialVec;  
% calculate the normal direction on each nodal point using the level set function %%%%%%%%%% 
tmpx = (Phi(:,[2:end 1])-Phi(:,[end 1:end-1]))/2/h;
tmpy = (Phi([2:end 1],:)-Phi([end 1:end-1],:))/2/h;
tmpx(:,1) = (Phi(:,2)-Phi(:,1))/h;   tmpx(:,end) = (Phi(:,end)-Phi(:,end-1))/h;    % Modify the boundary terms
tmpy(1,:) = (Phi(2,:)-Phi(1,:))/h;   tmpy(end,:) = (Phi(end,:)-Phi(end-1,:))/h;  
NormalVecX = -tmpx./sqrt(tmpx.^2+tmpy.^2 + 1e-16);                % add minus sign to make sure point outside 
NormalVecY = -tmpy./sqrt(tmpx.^2+tmpy.^2 + 1e-16);
NormalVec  = sqrt(tmpx.^2+tmpy.^2);
% Calculate <n*Sn>|_{\partial \Omega(t)}
nSn       = sigma11.*NormalVecX.^2 + 2*sigma12.*NormalVecX.*NormalVecY + sigma22.*NormalVecY.^2; 
perimeter = sum(sum(D_Phi.*NormalVec*h*h));                  % regularized length of the interface
nSn_ave   = sum(sum(D_Phi.*nSn.*NormalVec*h*h))/perimeter;   % regularized average of the driving force.
area_ratio = sum(sum(Hev(Phi)))/sum(sum(Hev(Phi0)));

delta_f   = (nSn - offset.*nSn_ave);    %CAUTION
Vm        = -beta*delta_f.^expN./(zz+zz1.*abs(delta_f).^expN) + lambda.*(1-area_ratio); 

%Vn        = Vm - gamma*RadialVec.*abs(erx.*NormalVecX + ery.*NormalVecY); 
%Vn         = Vm + gamma.*(erx0.*NormalVecX + ery0.*NormalVecY) - resist*(Xbart.*NormalVecX+Ybart.*NormalVecY); 
%Vn         =Vm + gamma.*(eee*erx0.*NormalVecX + ery0.*NormalVecY) + (gamma*resist/K)*(Xbart.*NormalVecX+Ybart.*NormalVecY)-0*NormalVecX;
Erx0 = (erx0+(eee-1)*(Xbart.*erx0+Ybart.*ery0).*Xbart);
Ery0 = (ery0+(eee-1)*(Xbart.*erx0+Ybart.*ery0).*Ybart);
%     Anisotropic velocity distribution stronger in direction of cell centroid velocity     
Vn         =Vm + gamma.*(Erx0.*NormalVecX + Ery0.*NormalVecY) + (gamma.*resist./K).*(Xbart.*NormalVecX+Ybart.*NormalVecY);
%     Anisotropic velocity distribution stronger in direction of cell centroid velocity    
%  The -1 in the last term accounts for lab frame
%resist=gamma*alpha/K 
 tp1 = Vm;
 tp2 = gamma.*(Erx0.*NormalVecX + Ery0.*NormalVecY) + (gamma.*resist./K).*(Xbart.*NormalVecX+Ybart.*NormalVecY);
 %tp2 = gamma.*(erx0.*NormalVecX + ery0.*NormalVecY) - (Xbart.*NormalVecX+Ybart.*NormalVecY);
% fprintf('AAAAAAAAAAAAAAAAA current area ratio is %.4f   \n',area_ratio);
