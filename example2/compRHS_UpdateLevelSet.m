% This code calculate the right-hand side of the level set function 
% \Phi_{t} = - Vn|\grad \Phi|;
% Zhiwen zhang/20130518

function [KPhi] = compRHS_UpdateLevelSet(t,h,Vn,Phi)
%% parameter
%% compute
% 1st & 2nd derivative terms
Dx_bj    = (Phi(:,1:end)     - Phi(:,[end 1:end-1]))/h;     % backward difference in x direction
Dx_fj    = (Phi(:,[2:end 1]) - Phi(:,1:end))/h;             % forward difference  in x direction
Dy_ib    = (Phi(1:end,:)     - Phi([end 1:end-1],:))/h;     % backward difference in y direction
Dy_if    = (Phi([2:end 1],:) - Phi(1:end,:))/h;             % forward difference  in y direction
DxmDxp_0 = (Phi(:,[2:end 1])-2*Phi(:,1:end)+Phi(:,[end 1:end-1]))/h^2;               % forward backward difference in x direction
DxmDxp_b = (Phi(:,1:end)-2*Phi(:,[end 1:end-1])+Phi(:,[end-1 end 1:end-2]))/h^2;       
DxmDxp_f = (Phi(:,[3:end 1 2])-2*Phi(:,[2:end 1])+Phi(:,1:end))/h^2;       
DymDyp_0 = (Phi([2:end 1],:)-2*Phi(1:end,:)+Phi([end 1:end-1],:))/h^2;               % forward backward difference in y direction
DymDyp_b = (Phi(1:end,:)-2*Phi([end 1:end-1],:)+Phi([end-1 end 1:end-2],:))/h^2;       
DymDyp_f = (Phi([3:end 1 2],:)-2*Phi([2:end 1],:)+Phi(1:end,:))/h^2;       
% Flux related terms
pxm = Dx_bj + 0.5*h*minmod(DxmDxp_0,DxmDxp_b);
pxp = Dx_fj - 0.5*h*minmod(DxmDxp_f,DxmDxp_0);
pym = Dy_ib + 0.5*h*minmod(DymDyp_0,DymDyp_b);
pyp = Dy_if - 0.5*h*minmod(DymDyp_f,DymDyp_0);

Pp  = sqrt(max(pxm,0).^2 + min(pxp,0).^2 + max(pym,0).^2 + min(pyp,0).^2);
Pm  = sqrt(min(pxm,0).^2 + max(pxp,0).^2 + min(pym,0).^2 + max(pyp,0).^2);
%% output:  
KPhi =  -(max(Vn,0).*Pp + min(Vn,0).*Pm);