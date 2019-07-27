% This code re-initialize the level set function into a sign-distance function
% Using same notation in Prof. Hou's 1999 JCP paper.  
% Using 2nd-ENO to discretize the Hamilton-Jacobi equation 
% CAUTION: Assume the level set function is periodic to deal with the Boundary Terms 
% So far only use the Forward Euler ODE solver; we may also consider the RK solver
% zhiwen zhang
% 2014/08/05
% Add diffusive term to test how this effects the result!!!
function  [Phi] = ReInitialization_LevelSetADDF(Phi,h,dt,eps_sign,eps_diff,Ntau) 
%% parameter
%% compute
sign_Phi = Phi./sqrt(Phi.^2+eps_sign.^2); 
% solve the re-initialize equation until steady-state
for i = 1:Ntau 
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
    % Update the level set function \Phi_{t} + sgn(\Phi^{*})|\grad \Phi| = sgn(\Phi^{*}) + eps_diff*Laplace(Phi); 
    Phi = Phi + dt*sign_Phi - dt*(max(sign_Phi,0).*Pp + min(sign_Phi,0).*Pm) + dt*eps_diff*(DxmDxp_0 + DymDyp_0);    
end
