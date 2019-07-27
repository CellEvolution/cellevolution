% This code generates the stiffness matrix. 
% Zhiwen Zhang/20170403
% Local Stencil pattern 
%  (7)(i-1,j+1)--- (3)(i,j+1)--- (6)(i+1,j+1)
%         |             |               |         
%  (4)(i-1,  j)--- (1)(i,  j)--- (2)(i+1,  j)
%         |             |               |
%  (8)(i-1,j-1)--- (5)(i,j-1)--- (9)(i+1,j-1)
% Input:  Nmuber of grids and parameters in the coefficient
% Output: Stiffness matrix and the load vector part (due to the Dirichlet BC)
% CAUTION: we have make the following mistakes. 
% (1) sign problems, (+ or -) in assign values for load vector b;
% (2) row_label: id used in the interior domain (free nodes); pt_label: id used in the fall domain (fixed nodes, call the BC values)
% CAUTION: I suspect that the above code is wrong, the index line 3-7 is wrong! (2014/07/31 find)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nx and Ny are the grids number along x- and y-axis
% CAUTION: in the following code, index i is for x corresponding to total Nx grids
%                                 index j is for y corresponding to total Ny grids 
function  [S,b] = RegularizedModel_StiffMat_Longtime(Nx,Ny,h,p,mu,lmd,BCx,BCy)
%% parameter %% Compute % Calculate the stiffness matrix 
% OLD CODEs
% spI = zeros(9*2*(N-2).^2,1);    % row    index
% spJ = zeros(9*2*(N-2).^2,1);    % column index
% spS = zeros(9*2*(N-2).^2,1);    % value  index
% b   = zeros(2*(N-2).^2,1);      % When the mode of the stencil touch BC, assign the BC values to the load vector b;
% NEW CODEs
spI = zeros(9*2*(Nx-2)*(Ny-2),1);    % row    index
spJ = zeros(9*2*(Nx-2)*(Ny-2),1);    % column index
spS = zeros(9*2*(Nx-2)*(Ny-2),1);    % value  index
b   = zeros(2*(Nx-2)*(Ny-2),1);      % When the mode of the stencil touch BC, assign the BC values to the load vector b;
Dummy = 2;                     % when spS=0, the dummy index has no effect, just to avoid negative values in coding.
%%%%%%%%%%% Interior Domain %%%%%%%%%%%%%%
for i = 3:Nx-2               
    for j = 3:Ny-2
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 = (i+1 -2)*(Ny-2)+j+1 -1;          
        tp7 = (i-1 -2)*(Ny-2)+j+1 -1;   tp8 = (i-1 -2)*(Ny-2)+j-1 -1;   tp9 = (i+1 -2)*(Ny-2)+j-1 -1;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    mu          lmd+2*mu   mu ...
                                                         (lmd+mu)/4   -(lmd+mu)/4  (lmd+mu)/4  -(lmd+mu)/4 ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          lmd+2*mu    mu         lmd+2*mu ...
                                                         (lmd+mu)/4   -(lmd+mu)/4  (lmd+mu)/4  -(lmd+mu)/4 ];
    end
end
%%%%%%%%%%% West BC %%%%%%%%%%%%%%
    for j = 3:Ny-2
        i = 2;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = Dummy;                    tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 = (i+1 -2)*(Ny-2)+j+1 -1;          
        tp7 = Dummy;                    tp8 =  Dummy;                   tp9 = (i+1 -2)*(Ny-2)+j-1 -1;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    mu          0          mu ...
                                                         (lmd+mu)/4    0           0           -(lmd+mu)/4 ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          lmd+2*mu    0          lmd+2*mu ...
                                                         (lmd+mu)/4    0           0           -(lmd+mu)/4 ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j;      %CAUTION: We need to get the label for current point                       
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc-h,ycc) +(lmd+mu)/4*BCy(xcc-h,ycc+h) - (lmd+mu)/4*BCy(xcc-h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc-h,ycc) +(lmd+mu)/4*BCx(xcc-h,ycc+h) - (lmd+mu)/4*BCx(xcc-h,ycc-h);
    end
%%%%%%%%%%% East BC %%%%%%%%%%%%%%
    for j = 3:Ny-2
        i = Nx-1;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 =  Dummy;                   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 =  Dummy;          
        tp7 = (i-1 -2)*(Ny-2)+j+1 -1;   tp8 = (i-1 -2)*(Ny-2)+j-1 -1;   tp9 =  Dummy;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)   0           mu          lmd+2*mu   mu ...
                                                          0           -(lmd+mu)/4  (lmd+mu)/4   0           ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)   0           lmd+2*mu    mu         lmd+2*mu ...
                                                          0           -(lmd+mu)/4  (lmd+mu)/4   0           ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc+h,ycc) -(lmd+mu)/4*BCy(xcc+h,ycc+h) + (lmd+mu)/4*BCy(xcc+h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc+h,ycc) -(lmd+mu)/4*BCx(xcc+h,ycc+h) + (lmd+mu)/4*BCx(xcc+h,ycc-h);
    end
%%%%%%%%%%% North BC %%%%%%%%%%%%%%
for i = 3:Nx-2               
        j = Ny-1;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 =  Dummy;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 =  Dummy;          
        tp7 =  Dummy;                   tp8 = (i-1 -2)*(Ny-2)+j-1 -1;   tp9 = (i+1 -2)*(Ny-2)+j-1 -1;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    0           lmd+2*mu   mu ...
                                                          0            0           (lmd+mu)/4  -(lmd+mu)/4 ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          0           mu         lmd+2*mu ...
                                                          0            0           (lmd+mu)/4  -(lmd+mu)/4 ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -       mu *BCx(xcc,ycc+h) -(lmd+mu)/4*BCy(xcc+h,ycc+h) + (lmd+mu)/4*BCy(xcc-h,ycc+h);
        b(2*row_label  ,1) = -(lmd+2*mu)*BCy(xcc,ycc+h) -(lmd+mu)/4*BCx(xcc+h,ycc+h) + (lmd+mu)/4*BCx(xcc-h,ycc+h);
end
%%%%%%%%%%% South BC %%%%%%%%%%%%%%
for i = 3:Nx-2               
        j = 2;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 =  Dummy;                   tp6 = (i+1 -2)*(Ny-2)+j+1 -1;          
        tp7 = (i-1 -2)*(Ny-2)+j+1 -1;   tp8 =  Dummy;                   tp9 =  Dummy;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    mu          lmd+2*mu   0 ...
                                                         (lmd+mu)/4   -(lmd+mu)/4  0           0           ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          lmd+2*mu    mu         0 ...
                                                         (lmd+mu)/4   -(lmd+mu)/4  0           0           ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -       mu *BCx(xcc,ycc-h) + (lmd+mu)/4*BCy(xcc+h,ycc-h) - (lmd+mu)/4*BCy(xcc-h,ycc-h);
        b(2*row_label  ,1) = -(lmd+2*mu)*BCy(xcc,ycc-h) + (lmd+mu)/4*BCx(xcc+h,ycc-h) - (lmd+mu)/4*BCx(xcc-h,ycc-h);
end
%%%%%%%%%%% North East Corner %%%%%%%%%%%%%%
        i = Nx-1; j = Ny-1;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 =  Dummy;                   tp3 = Dummy;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 = Dummy;          
        tp7 =  Dummy;                   tp8 = (i-1 -2)*(Ny-2)+j-1 -1;   tp9 = Dummy;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  0            0          lmd+2*mu   mu ...
                                                          0            0           (lmd+mu)/4  0           ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  0            0          mu         lmd+2*mu ...
                                                          0            0           (lmd+mu)/4  0           ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc+h,ycc)   -        mu *BCx(xcc,ycc+h) ...
                             -(lmd+mu)/4*BCy(xcc+h,ycc+h) + (lmd+mu)/4*BCy(xcc-h,ycc+h) + (lmd+mu)/4*BCy(xcc+h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc+h,ycc)   - (lmd+2*mu)*BCy(xcc,ycc+h) ...
                             -(lmd+mu)/4*BCx(xcc+h,ycc+h) + (lmd+mu)/4*BCx(xcc-h,ycc+h) + (lmd+mu)/4*BCx(xcc+h,ycc-h);
%%%%%%%%%%% North West Corner %%%%%%%%%%%%%%
        i = 2; j = Ny-1;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 = Dummy;          
        tp4 = Dummy;                    tp5 = (i   -2)*(Ny-2)+j-1 -1;   tp6 = Dummy;          
        tp7 = Dummy;                    tp8 = Dummy;                    tp9 = (i+1 -2)*(Ny-2)+j-1 -1;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    0          0         mu ...
                                                          0            0           0        -(lmd+mu)/4  ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          0          0         lmd+2*mu ...
                                                          0            0           0        -(lmd+mu)/4  ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc-h,ycc)   -        mu *BCx(xcc,ycc+h) ...
                             -(lmd+mu)/4*BCy(xcc+h,ycc+h) + (lmd+mu)/4*BCy(xcc-h,ycc+h) - (lmd+mu)/4*BCy(xcc-h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc-h,ycc)   - (lmd+2*mu)*BCy(xcc,ycc+h) ...
                             -(lmd+mu)/4*BCx(xcc+h,ycc+h) + (lmd+mu)/4*BCx(xcc-h,ycc+h) - (lmd+mu)/4*BCx(xcc-h,ycc-h);
%%%%%%%%%%% South West Corner %%%%%%%%%%%%%%
        i = 2; j = 2;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = (i+1 -2)*(Ny-2)+j   -1;   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = Dummy;                    tp5 = Dummy;                    tp6 = (i+1 -2)*(Ny-2)+j+1 -1;          
        tp7 = Dummy;                    tp8 = Dummy;                    tp9 = Dummy;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  lmd+2*mu    mu          0     0 ...
                                                         (lmd+mu)/4    0           0           0     ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  mu          lmd+2*mu    0     0 ...
                                                         (lmd+mu)/4    0           0           0     ];
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc-h,ycc)   -        mu *BCx(xcc,ycc-h) ...
                             +(lmd+mu)/4*BCy(xcc-h,ycc+h) - (lmd+mu)/4*BCy(xcc-h,ycc-h) + (lmd+mu)/4*BCy(xcc+h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc-h,ycc)   - (lmd+2*mu)*BCy(xcc,ycc-h) ...
                             +(lmd+mu)/4*BCx(xcc-h,ycc+h) - (lmd+mu)/4*BCx(xcc-h,ycc-h) + (lmd+mu)/4*BCx(xcc+h,ycc-h);
%%%%%%%%%%% South East Corner %%%%%%%%%%%%%%
        i = Nx-1; j = 2;
        row_label = (i-2)*(Ny-2)+j-1;     % transform the global index to linear equation index
        spI((row_label-1)*18 +1:(row_label-1)*18 +9,1) = 2*row_label-1;    % row   label for first  Eq, laplace{u} 
        spI((row_label-1)*18+10:(row_label-1)*18+18,1) = 2*row_label;      % row   label for second Eq, laplace{v}
        % column label   
        tp1 = (i   -2)*(Ny-2)+j   -1;   tp2 = Dummy;   tp3 = (i   -2)*(Ny-2)+j+1 -1;          
        tp4 = (i-1 -2)*(Ny-2)+j   -1;   tp5 = Dummy;   tp6 = Dummy;          
        tp7 = (i-1 -2)*(Ny-2)+j+1 -1;   tp8 = Dummy;   tp9 = Dummy;          
        spJ((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[2*tp1-1 2*tp2-1 2*tp3-1 2*tp4-1 2*tp5-1 2*tp6   2*tp7   2*tp8   2*tp9];   %odd for u, even for v;                    
        spJ((row_label-1)*18+10:(row_label-1)*18+18,1)=[2*tp1   2*tp2   2*tp3   2*tp4   2*tp5   2*tp6-1 2*tp7-1 2*tp8-1 2*tp9-1]; %odd for u, even for v;                    
        % matrix value label
        spS((row_label-1)*18 +1:(row_label-1)*18 +9,1)=[-(2*lmd+6*mu)  0            mu          lmd+2*mu   0 ...
                                                          0           -(lmd+mu)/4   0           0          ];
        spS((row_label-1)*18+10:(row_label-1)*18+18,1)=[-(2*lmd+6*mu)  0            lmd+2*mu    mu         0 ...
                                                          0           -(lmd+mu)/4   0           0          ];                                                   
        % assign the corresponding values to the load vector
        pt_label = i*Ny + j; 
        xcc = p(pt_label,1);  ycc = p(pt_label,2);     % coordinates of the current center
        b(2*row_label-1,1) = -(lmd+2*mu)*BCx(xcc+h,ycc)   -        mu *BCx(xcc,ycc-h) ...
                             -(lmd+mu)/4*BCy(xcc+h,ycc+h) - (lmd+mu)/4*BCy(xcc-h,ycc-h) + (lmd+mu)/4*BCy(xcc+h,ycc-h);
        b(2*row_label  ,1) = -       mu *BCy(xcc+h,ycc)   - (lmd+2*mu)*BCy(xcc,ycc-h) ...
                             -(lmd+mu)/4*BCx(xcc+h,ycc+h) - (lmd+mu)/4*BCx(xcc-h,ycc-h) + (lmd+mu)/4*BCx(xcc+h,ycc-h);
%% output
S = sparse(spI,spJ,spS,2*(Nx-2)*(Ny-2),2*(Nx-2)*(Ny-2));
