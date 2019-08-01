% This code use the finite difference method and level set method to solve
% the cell evolution model of Zhang et al.
% EARLIER COMMENTS:
% The elastic model is linear elasticity PLANE-STRAIN
% Zhiwen Zhang / 20141124
% Convergence study and area preservation study!
% On the Interface, replace Heaviside function by 1/2 to avoid numerical difficulty
% Local Stencil pattern 
%  (7)(i-1,j+1)--- (3)(i,j+1)--- (6)(i+1,j+1)
%         |             |        t       |         
%  (4)(i-1,  j)--- (1)(i,  j)--- (2)(i+1,  j)
%         |             |               |
%  (8)(i-1,j-1)--- (5)(i,j-1)--- (9)(i+1,j-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Choose Lx > Ly so that we can observe longtime evolution along x-axis.
%  Therefore, all the codes related to L need to change. 
%  Zhiwen Zhang-----------2017-4-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Recompute the solution using finer mesh
%  Zhiwen Zhang-----------2018-07-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% parameter
mu = 1; lmd = 1;                    % mu/lmd: lame constant;  This setting is equivalent to (poisson ratio) nu = 1/4; (Young's modulus) E = 1;
everystep = 1;                      % reinitialize every ? step 
%sigma  = 0.0;                       % always zero for BC
beta   =2.5;                        % control the amplitude of the nonlinear velocity in G(z)
eee = 3;                            %ellipticity coefficient in velocity distribution; eee=1 is isotropic circular symmetry
K  = 3;                            % amplitute used in computing the self-stress of Cell. 
lambda = -30;
gamma  = 0.8;
%a = 1;                               % the radius of circular cells
h = 1/16;                             % mesh size in physical domain partition
%% %%%%%%%%%%%%% parameters used in regularization and reinitialization of level set function
r    = 3.5;                              % regularization ratio 
eps      = r*h;                          % parameter used in the regularized delta function 
eps_sign = 1*h;                          % parameter used in the smoothed sigh function
eps_diff = h*h*h/2;                     % parameter used in the diffusive term when reinitialize level set! ?
%% %%%%%%%%%%%%%% parameters used in setting mesh grid  and time-step %%%%%%%%%%%%%%
filename  = 'DataDuroturning_Ly10h1over16_tf32_eee3_resist0';        % adaptively compute self-stress tensor
filename2 = 'DataDuroturning_Ly10h1over16_tf32_eee3_resist0_ext';    % adaptively compute self-stress tensor
figure_id0 = 217;
% snapshot time data 
dt = h/20; t0 = 0.00; tf = 32;  tt = t0:dt:tf-dt;  nnt = length(tt);     % computational time, sigma = +2.0; 
tss=320;  itv=round(nnt/tss);  
% generate the fine mesh grid and related data structure  
Lx = 16;  Ly = 10;  [X,Y] = meshgrid(-Lx:h:Lx,-Ly:h:Ly);       % NEW CODE for rectangular domain
p           = [X(:),Y(:)];     
Nx = length(-Lx:h:Lx);  Ny = length(-Ly:h:Ly);
Phisnap=zeros(Ny,Nx,tss);  s=zeros(nnt,1);                       
nSnsnap=zeros(Ny,Nx,tss);   
Usnap  =zeros(Ny,Nx,tss);  Vsnap=zeros(Ny,Nx,tss);   
xbars = zeros(tss,1);      ybars = zeros(tss,1);
dragXsnap = zeros(Ny,Nx,tss);    dragYsnap = zeros(Ny,Nx,tss); 
%% parameters used in label the grid in FEM and BC %%%%%%%%%%%%%%%%%%%%%%%%
p0          = [0,0];
NP          = size(p,1);                               % All nodes
idx1        = find(abs(p(:,1)-p0(1))>Lx-0.1*h | abs(p(:,2)-p0(2))>Ly-0.1*h);
idx2        = find(abs(p(:,1)-p0(1))<Lx-0.1*h & abs(p(:,2)-p0(2))<Ly-0.1*h);     % check and label fine grid BC
[Xint,Yint] = meshgrid(-Lx+h:h:Lx-h,-Ly+h:h:Ly-h);  
pint        = [Xint(:),Yint(:)];   
NPint       = size(pint,1);                            % Free nodes (exclude the BC nodes)
rho         = sqrt(X.^2 + Y.^2 + 1e-16);               % polar distance 
% generate the Dirichlet BC function 
% (2) BC used for cell evolution!!!
DBCx        = @(x,y) 0;
DBCy        = @(x,y) 0;
% generate the level set function and the sign function according to the initial cell interface 
tol            = 0.001; 
Ntau           = 40;
% generate the regularized delta and heviside function 
Delta_eps       = @(z) (abs(z)<eps)*0.5.*(1+cos(pi*z/eps))/eps; 
Hev_eps         = @(z) (abs(z)<=eps).*((z+eps)/2/eps + sin(pi*z/eps)/2/pi) + (z>eps)*1; 
AvgHev_eps      = @(z)  1/2;                                  % averaged heviside function   
%  initial cell shape parameters
X0 = -5;  
a0 = 1.0;
b0 = 3.4;
n0 = 2; 
r0 = 1.8;
c0=1; c1=1;
rr1  = ((sqrt((X/c0-X0).^2+((Y-3)/c1).^2)-r0)/a0).^n0;
rr2  = (atan2((Y-3)/c1,X/c0-X0)/b0).^n0;
Phi  = -1*Hev_eps(6-abs(atan2((Y-3)/c1,X/c0-X0))).*(rr1 + rr2 - 1)-0.001;
Phi0 = Phi;   % record the initial shape 
% Phisnap(:,:,1) = Phi0; 
    vv = [0,0];
    H_Phi     = Hev_eps(Phi);            D_Phi     = Delta_eps(Phi);
    Xbar0 = sum(sum(Hev_eps(Phi).*X))/sum(sum(Hev_eps(Phi)));   % centroid of the cell
    Ybar0 = sum(sum(Hev_eps(Phi).*Y))/sum(sum(Hev_eps(Phi)));

%% Auxiliary Variables for plot and debug
id_i = 1;      vali(id_i,1) = 1;                   vali(id_i,2) = 1;                        % time
id_i = 2;      vali(id_i,1) = vali(id_i-1,2)+1;    vali(id_i,2) = vali(id_i,1) + 0;         % area ratio
vals = zeros(10, vali(end,end));
id_j    = 1;
%% compute 
% Calculate the stiffness matrix and load vector (minus body force). CAUTION: this part is INVARIANT with respect to time
[S,b0] = RegularizedModel_StiffMat_Longtime(Nx,Ny,h,p,mu,lmd,DBCx,DBCy); 
Fixed  = [2*idx1-1; 2*idx1];   Free = sort([2*idx2-1; 2*idx2]);    % use 'sort' to comply with the order of unknowns in the stiffness matrix
% begin the time evolution  CAUTION: The code is ReInitialization_LevelSetADDF.m is reusable. 
[Phi] = ReInitialization_LevelSetADDF(Phi,h,dt,eps_sign,eps_diff,10*Ntau);   % ADd DiFfusive term; 
snapID = 0;  
Xbart = 0;         Ybart = 0.0;      
for n = 1:nnt 
    t = tt(n); 
    b         = b0;    
    % calculate the variant part of the load vector (own to the decay of the delta function, we only consider the interior domain)
    % Compute Extra Term 
    Phix=(Phi(:,[2:end 1])-Phi(:,[end 1:end-1]))/2/h;  Phiy=(Phi([2:end 1],:)-Phi([end 1:end-1],:))/2/h;
    Phix(:,1) = (-Phi(:,3)+4*Phi(:,2)-3*Phi(:,1))/2/h;   Phix(:,end) = (3*Phi(:,end)-4*Phi(:,end-1)+Phi(:,end-2))/2/h;    % Modify the boundary terms
    Phiy(1,:) = (-Phi(3,:)+4*Phi(2,:)-3*Phi(1,:))/2/h;   Phiy(end,:) = (3*Phi(end,:)-4*Phi(end-1,:)+Phi(end-2,:))/2/h; 
    H_Phi     = Hev_eps(Phi);            D_Phi     = Delta_eps(Phi);
    Xbar = sum(sum(Hev_eps(Phi).*X))/sum(sum(Hev_eps(Phi)));   % centroid of the cell
    Ybar = sum(sum(Hev_eps(Phi).*Y))/sum(sum(Hev_eps(Phi)));
    if n==1 
       Xbar0 = Xbar;                    Ybar0 = Ybar;         
    end
    Xbart = (Xbar-Xbar0)/dt;         Ybart = (Ybar-Ybar0)/dt;     
    Xbar0 = Xbar;                    Ybar0 = Ybar; 
    Zx   = X-Xbar;          Zy   = Y-Ybar; 
     tmp1 = (-H_Phi).*(K*(Zx+(eee-1).*(Xbart.*Zx+Ybart.*Zy).*Xbart)); 
     tmp2 = (-H_Phi).*(K*(Zy+(eee-1).*(Xbart.*Zx+Ybart.*Zy).*Ybart));
%     Anisotropic velocity distribution stronger in direction of cell centroid velocity     
    for i = 3:Nx-2         % i is colume index; j is row index;      
        for j = 3:Ny-2
            row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
            % assign the corresponding values to the load vector                   
            b(2*row_label-1,1) = b(2*row_label-1,1) - h*h*tmp1(j,i);  % 
            b(2*row_label  ,1) = b(2*row_label  ,1) - h*h*tmp2(j,i);  %              
        end
    end
    % solve the linear equation to obtain (u,v)
    UV = zeros(2*Nx*Ny,1);       % displacement solution [u1 v1 u2 v2,...]';
    UV(2*idx1-1,1)  = DBCx(p(idx1,1),p(idx1,2)); 
    UV(2*idx1   ,1) = DBCy(p(idx1,1),p(idx1,2));     
    UV(Free)        = S\b;
    u1D = UV(1:2:end-1,1);   v1D = UV(2:2:end,1); 
    U   = reshape(u1D,Ny,Nx);  V   = reshape(v1D,Ny,Nx); 
    % calculate the stresss tensor      
      
      [Vn,NormalVecX,NormalVecY,nSn,nSn_ave,Ux,Uy,Vx,Vy,tp1,tp2,sigma11,sigma22,sigma12] = ...
          NormalVelocity(X,Y,Phi0,Phi,U,V,p,mu,lmd,beta,eee,lambda,K,gamma,Hev_eps,Delta_eps,Xbart,Ybart);   
  
     % CAUTION: the above kinetic code for velocity is reusable!!!
    [KPhi1] = compRHS_UpdateLevelSet(t,       h,Vn,Phi             );
    [KPhi2] = compRHS_UpdateLevelSet(t+0.5*dt,h,Vn,Phi+0.5*dt*KPhi1);
    [KPhi3] = compRHS_UpdateLevelSet(t+0.5*dt,h,Vn,Phi+0.5*dt*KPhi2);
    [KPhi4] = compRHS_UpdateLevelSet(t+1.0*dt,h,Vn,Phi+1.0*dt*KPhi3);
     Phi    = Phi  + (dt/6) * (KPhi1  + 2*KPhi2  + 2*KPhi3  + KPhi4); 
     t      = t + dt;      
    %% ReInitialization the Level Set function
    if mod(n,everystep)==0 
      [Phi] = ReInitialization_LevelSetADDF(Phi,h,dt,eps_sign,eps_diff,Ntau);   % ADd DiFfusive term; 
      fprintf('At time t=%.4f, do reInitialization once. \n',t);
    end
    %% compute area and ratio 
    s(n,1)              = sum(sum(Hev_eps(Phi)))*h*h; 
    s0                  = sum(sum(Hev_eps(Phi0)))*h*h;     
    % time
    vals(id_j,1) = t;
    vals(id_j,vali(2,1):vali(2,2)) = s(n,1)/s0; 
    %% save the snapshot solution 
    if mod(n,itv) == 0;
       snapID                = snapID + 1;
       Phisnap(:,:,snapID)   = Phi; 
       dragXsnap(:,:,snapID) = tmp1;
       dragYsnap(:,:,snapID) = tmp2;
       nSnsnap(:,:,snapID) = nSn;
       Usnap(:,:,snapID)   = U;     Vsnap(:,:,snapID)   = V;
       xbars(snapID,1)     = Xbar;  ybars(snapID,1)     = Ybar; 
       fprintf(1,'Nx=%4d, Ny=%4d, Ntau=%2d, finish update time t=%.4f of total tf=%.4f \n',Nx,Ny,Ntau,tt(n),tf);
       save(filename,'tt','s','Phisnap','dragXsnap','dragYsnap','X','Y','Lx','Ly','Phi0','xbars','ybars');    
       save(filename2,'tt','s','nSnsnap','Usnap','Vsnap','X','Y','Lx','Ly','Phi0','xbars','ybars');    

    end

    id_j = id_j + 1;
end
