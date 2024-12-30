function  PreDE
% Preconditioning of 2D Heat Equation with Convection
% Lagrange FEM, homogeneous Dirichlet boundary condition
% domain (0,2pi) x (0,2pi)
% - preconditioning
% - steps of CG and GMRES
% - spectral bounds and their graphs

Ntri = 24 % number of intervals in both x1 and x2
Convection = 1 % if convection term (1) or not (0) (then GMRES or CG)
Time = 1 % time evolution (1) or not (0)
ShowSpec = 1 % plot spectrum (1) or not (0)
if (Convection==1 && Time==0) disp('Sorry, Time must be 1 if Convection==1.'); return; end;

E = eye(2); % for homogenization
toler = 1e-9; % tolerance for all methods
max_steps = 1000; % max steps for all methods
L1 = 2*pi; L2 = 2*pi; % edges of the domain domain 
NNN = ones(2,1)*Ntri; % number of intervals 
N1 = NNN(1); N2 = NNN(2); % number of intervals 
[N1,N2]
h1 = L1/N1; h2 = L2/N2; 
dt = 0.5 % time step
Nnod = (N1+1)*(N2+1); % number of all nodes
Nvox = N1*N2;  % number of pixels
Nele = 2*Nvox;  % number of elements
Nnod_in = (N1-1)*(N2-1); % number of all inner nodes

[nodes,peri,marg] = Nodes(N1,N2,L1,L2); % all nodes including whole boundary

elem = Elem(N1,N2); % triangles of FEM
A = zeros(1*Nnod); % stiffness matrix
AC = zeros(Nnod); % for convectine term
PR = zeros(Nnod); % preconditioning matrix
PZS = zeros(Nnod); % preconditioning matrix
AT = zeros(Nnod); % mass matrix
B = zeros(1*Nnod,1); % rhs
D1 = zeros(Nele,Nnod); % 2 matrices of derivatives in quadrature points (dx and dy)
D2 = zeros(Nele,Nnod);

alfa1 = zeros(Nnod,1)+1e12;
alfa2 = zeros(Nnod,1)-1;
beta = 0;

Fpom_a = Fxy_data_a(1,1); % class for conductivity a(x,y) 
Fpom_f = Fxy_RHS(1,1); % class for RHS
Fpom_c = Fxy_data_c(1,1); % class for convection c(x,y) 
for kk = 1:Nele
    kde = elem(kk,:);
    nod = nodes(kde,:);  
    pom = sum(nod)/3;
    s1 = pom(1);
    s2 = pom(2);
    dd = abs(det([nod(1,:)-nod(3,:);nod(2,:)-nod(3,:)]))/2;
    nod1 = [nod,ones(3,1)];
    der = inv(nod1);
    der = der(1:2,:); % d/dx, d/dy for 3 basis functions
    g = sum(nod)/3; % mass center    
    a = Fpom_a.Eval(s1,s2); 
   
    apre = Fpom_a.Fxy_ref_a;
    aq = inv(a);    
    fc = Fpom_c.Eval(s1,s2); 
    pomb = ones(3,1)*fc*der*dd/3;   % CONVECTIVE TERM
    f = Fpom_f.Eval(s1,s2);    
    pom8 = [kk,kk+1*Nele];
    CA(pom8,pom8) = a;      
    CAq(pom8,pom8) = aq;
    
    % derivatives [dxu1,0,dyu1]' and [0,dyu2,dxu2]' in all quad.points:
    derS = [der(1,:),zeros(1,3);
        zeros(1,3),der(2,:);
        der(2,:),der(1,:)];  
    derr = [der(1,:);der(2,:)]; 
    pom = derr'*a*derr*dd;
    pompre = derr'*apre*derr*dd ;
    pomtime = 1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/dt;

    D1(kk,kde) = der(1,:); 
    D2(kk,kde) = der(2,:);

    A([kde],[kde]) = A([kde],[kde]) + pom ;
    AC([kde],[kde]) = AC([kde],[kde]) + pomb;     
    AT([kde],[kde]) = AT([kde],[kde]) + pomtime;
    PR([kde],[kde]) = PR([kde],[kde]) + derr'*derr*dd ; % Green pre
    PZS([kde],[kde]) = PZS([kde],[kde]) + derr'*apre*derr*dd ; % Green anizo pre
    B([kde]) = B([kde])+f*dd/3;   

    pomb2 = (pomb'-pomb)/2; 
    p2 = max(imag(eig(pinv(pompre+pomtime*Time)*(pomb2))));
    p = sort(real(eig(pinv(pompre+pomtime*Time)*(pom+pomtime*Time))));
    if (beta<p2) beta = p2; end;
    if (Time==1) kkk = 1; else kkk = 2; end;
    for j = 1:3
        s = kde(j);        
        if (p(kkk)<alfa1(s)) alfa1(s) = p(kkk); end;
        if (p(3)>alfa2(s)) alfa2(s) = p(3); end;
    end
end;

% AT0 = AT;
% PR0 = PR;

if (Time==1)
D_AP = diag(PR+AT*Time)./diag(A+AT*Time);
D_APZS = diag(PZS+AT*Time)./diag(A+AT*Time);

% arrays for lower and upper bounds to eienvalues of preconditioned matrix:
alfa1_Jac = zeros(Nnod,1)+1e9;
alfa2_Jac = zeros(Nnod,1)-1;
alfa1_JacZS = zeros(Nnod,1)+1e9;
alfa2_JacZS = zeros(Nnod,1)-1;
betaD = 0;
for kk = 1:Nele
    kde = elem(kk,:);
    nod = nodes(kde,:);  
    pom = sum(nod)/3;
    s1 = pom(1);
    s2 = pom(2);
    dd = abs(det([nod(1,:)-nod(3,:);nod(2,:)-nod(3,:)]))/2;
    nod1 = [nod,ones(3,1)];
    der = inv(nod1);
    der = der(1:2,:); % d/dx, d/dy for 3 basis functions
    g = sum(nod)/3; % mass center    
    a = Fpom_a.Eval(s1,s2); 
    apre = Fpom_a.Fxy_ref_a;
    % aq = inv(a);   
    pom8 = [kk,kk+1*Nele];     
   
    derr = [der(1,:);der(2,:)]; 
    pom = derr'*a*derr*dd;
    pompre = derr'*apre*derr*dd ;
    pomtime = 1*[2,1,1;1,2,1;1,1,2]*h1*h2/24/dt;

    pomb2 = (pomb'-pomb)/2; 
    p2 = max(imag(eig(pinv(pompre+pomtime)*(pomb2+pomtime))));    

    pomdiag2 = diag(real(sqrt(D_AP(kde))));
    pomdiag2ZS = diag(real(sqrt(D_APZS(kde))));

    pB = sort(real(eig(pinv(derr'*derr*dd+pomtime*Time)*(pomdiag2*(pom+pomtime*Time)*pomdiag2))));
    pBZS = sort(real(eig(pinv(pompre+pomtime*Time)*(pomdiag2ZS*(pom+pomtime*Time)*pomdiag2ZS))));
    
    pAC_BZS = sort(imag(eig(pinv(pompre+pomtime*Time)*(pomdiag2ZS*(pomb2)*pomdiag2ZS))));
    p2 = max(pAC_BZS);
    if (betaD<p2) betaD = p2; end;
    if (Time==1) kkk = 1; else kkk = 2; end;  % only for Time==1 
    for j = 1:3
        s = kde(j);        
        if (pB(kkk)<alfa1_Jac(s)) alfa1_Jac(s) = pB(kkk); end
        if (pB(3)>alfa2_Jac(s)) alfa2_Jac(s) = pB(3); end
        if (pBZS(kkk)<alfa1_JacZS(s)) alfa1_JacZS(s) = pBZS(kkk); end
        if (pBZS(3)>alfa2_JacZS(s)) alfa2_JacZS(s) = pBZS(3); end
    end
end;
end % of computing bounds for Time==1

AR = A;
BR = B;

D = [D1;D2];
PR = D'*D*dd;
poma = Fpom_a.Fxy_ref_a;
% PZS = [poma(1,1)*D1'*D1+poma(2,2)*D2'*D2+poma(1,2)*D1'*D2+poma(2,1)*D2'*D1]*dd;

% applying boundary conditions and removing some bounds:
for j = Nnod:-1:1
    if (marg(j)>0) 
        AR(:,j) = [];
        AR(j,:) = [];
        AC(:,j) = [];
        AC(j,:) = [];
        AT(:,j) = [];
        AT(j,:) = [];
        BR(j) = [];        
        PR(:,j) = [];
        PR(j,:) = [];
        PZS(:,j) = [];
        PZS(j,:) = [];
        alfa1(j) = [];
        alfa2(j) = [];
        if (Time==1)            
        alfa1_Jac(j) = [];
        alfa2_Jac(j) = [];
        alfa1_JacZS(j) = [];
        alfa2_JacZS(j) = [];
        end;
    end
end


% CG:

if (Convection==0)
[U1,stFEM1] = CGP_3_without(AR+AT*Time,BR,BR*0,toler,max_steps);
[x01,stFEM_pre1,ppp,toc1] = CGP_3_pre(AR+AT*Time,BR,BR*0,toler,max_steps,pinv(PZS+AT*Time));
stepsCG_ZS_Axb_n_p = [stFEM1,stFEM_pre1]

pom = (diag(AR+AT*Time)./diag(PR+AT*Time)).^(1/2);
Da = diag(pom);
PPD = Da*(PR+AT*Time)*Da;
[x01,stFEM_pre1,ppp,toc1] = CGP_3_pre(AR+AT*Time,BR,BR*0,toler,max_steps,pinv(PPD));
stepsCG_Berta_Axb_n_p = [stFEM_pre1]
Ux = U1;

STEPS_CG = [stepsCG_ZS_Axb_n_p,stepsCG_Berta_Axb_n_p]%,stepsCG_ZSBerta_Axb_n_p]

end % of Convection==0



% GMRES:

if (Convection==1)
ARC = AR+AC;

[U1,f,r,stFEM1] = gmres(ARC+AT*Time,BR,[],toler,max_steps);
[x01,f,r,stFEM_pre1] = gmres(ARC+AT*Time,BR,[],toler,max_steps,PZS+AT*Time);
stepsGM_ZS_Axb_n_p = [stFEM1,stFEM_pre1]
Ux = U1;

pom = (diag(ARC+AT*Time)./diag(PR+AT*Time)).^(1/2);
Da = diag(pom);
PPD = Da*(PR+AT*Time)*Da;
[U1,f,r,stFEM1] = gmres(ARC+AT*Time,BR,[],toler,max_steps);
[x01,f,r,stFEM_pre1] = gmres(ARC+AT*Time,BR,[],toler,max_steps,PPD);
stepsGM_Berta_Axb_n_p = [stFEM_pre1]
Ux = U1;

STEPS_GMRES = [stepsGM_ZS_Axb_n_p,stepsGM_Berta_Axb_n_p]%,stepsGM_ZSBerta_Axb_n_p]

end % of Convection==1


% plots of spectra and bounds:

if (ShowSpec==1 && Convection==0)    

subplot(1,2,1);
p = sort(real(eig(pinv(PZS+AT*Time)*(AR+AT*Time))));
cla; hold on;
plot(sort(alfa1),'r.-')
plot(sort(alfa2),'g.-')
title('Green pre')
if (Time==0) axis([0,Nnod_in,0,max(alfa2)])
else axis([0,Nnod_in,0,max(max(alfa2_Jac),max(alfa2))])
end;
plot(p,'b.')
ConditionAZS = max(p)/min(p);
pom11 = max(alfa2);

subplot(1,2,2);
cla; hold on;
pom = (diag(AR+AT*Time)./diag(PR+AT*Time)).^(1/2);
Da = diag(1./pom);
PPD = PR+AT*Time;
p22 = sort(real(eig(pinv(PPD)*(Da*(AR+AT*Time)*Da))));
if (Time==1)
plot(sort(alfa1_Jac),'r.-')
plot(sort(alfa2_Jac),'g.-')
end;
title('Green+Jacobi pre')
if (Time==0) axis([0,Nnod_in,0,max(alfa2)])
else axis([0,Nnod_in,0,max(max(alfa2_Jac),max(alfa2))])
end;
ConditionABer = max(p22)/min(p22); 
plot(p22,'b.')

end %  of Convection==0
    







if (ShowSpec==1 && Convection==1) % Convection==1

subplot(1,2,1);
cla; hold on;
p777 = eig(pinv(PZS+AT*Time)*(ARC+AT*Time));
m1 = min(alfa1);
m2 = max(alfa2);
p = beta;
plot([m1,m2,m2,m1,m1],[-p,-p,p,p,-p],'r-','LineWidth',1);
plot(real(p777),imag(p777),'b.')
title('Green pre')

subplot(1,2,2);
cla; hold on;
pom = (diag(ARC+AT*Time)./diag(PR+AT*Time)).^(1/2);
Da = diag(1./pom);
PPDa = Da*(ARC+AT*Time)*Da;
p = eig(pinv(PR+AT*Time)*PPDa);
plot(real(p),imag(p),'b.')
title('Green+Jacobi pre')

if (Time==1)
    plot([min(alfa1_JacZS),max(alfa2_JacZS),max(alfa2_JacZS),min(alfa1_JacZS),min(alfa1_JacZS)],[-betaD,-betaD,betaD,betaD,-betaD],'r-','LineWidth',1)

    z1 = max([max(alfa2),max(alfa2_JacZS)]);
    z2 = max([betaD,beta])*1.2;
    subplot(1,2,1)
    axis([0,z1,-z2,z2]);
    subplot(1,2,2)
    axis([0,z1,-z2,z2]);      
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function elem = Elem(N1,N2)
NN = N1*N2;
elem = zeros(2*NN,3);
j = 1;
for k2 = 1:N2
    for k1 = 1:N1
        p = (k2-1)*(N1+1)+k1; % lower left corner
            elem(j,:) = [p,p+1,p+N1+2]; 
            elem(j+1,:) = [p,p+N1+1,p+N1+2];
            j = j+2;
    end
end

function [nodes,peri,marg] = Nodes(N1,N2,L1,L2)
NN = (N1+1)*(N2+1);
h1 = L1/N1; h2 = L2/N2; 
nodes = zeros(NN,2);
peri = zeros(NN,1);
orig = zeros(NN,1);
j = 1;
for k2 = 1:N2+1
    for k1 = 1:N1+1        
            nodes(j,:) = [(k1-1)*h1,(k2-1)*h2];
            if (k1>N1+0.5) peri(j) = 1; end
            if (k2>N2+0.5) peri(j) = 2; end
            if (k1>N1+0.5 && k2>N2+0.5 ) peri(j) = 3; end
            if (peri(j)>0 || k1==1 || k2==1) marg(j) = 1; end
            j = j+1;
    end
end

function y = FindNeighElem(elem,N)
y = zeros(N,6); % finds neighboring elements
z = zeros(N,1);
for k = 1:N
    for r = 1:size(elem,1)
    for s = 1:3
        if (elem(r,s)==k)
            z(k) = z(k) + 1 ;
            y(k,z(k)) = r;
        end
    end
    end
end

function u = FindNeighNodes(elem,N)
y = FindNeighElem(elem,N);
u = zeros(N,9); % neighboring nodes
for k = 1:N % node number k
    pom = 0;
    for s = 1:6 % elements attached to k-th node       
        if (y(k,s)>0)
            pom = [pom,elem(y(k,s),:)];
        end  
    end
    pom = sort(pom);
    for j = size(pom,2):-1:2
        if (pom(j)==pom(j-1) || pom(j)==k)
            pom(j)=-1;
        end            
    end   
    pom = sort(pom);
    pom = pom(find(pom>0)) ;   
    u(k,1:size(pom,2)) = pom;    
end





