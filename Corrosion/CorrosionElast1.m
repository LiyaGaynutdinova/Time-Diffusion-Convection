% FEM for 2d problem
% -div.(a grad u)=f on (0,L1)x(0,L2), 
% periodic BC
function abc
clc;

L1 = 1; % length of domain in x1
L2 = 1; % length of domain in x2
N1 = 8; % EVEN number of intervals in x1
N2 = 8; % EVEN number of intervals in x2
h1 = L1/N1;
h2 = L2/N2;
max_steps = 3000;
gc_mac = [1e-6;1e-6]; % macro gradient
gc_mac = [1;1]*1e-5;
eps_mac = [1;1;0]*1e-5;
c_mac = 2e-5;
threshold_c = 1e-9; % 1e-9
molar_mass = 89;  % 89
density = 0.8e-3; % 3.4e-3
Gc = 1e-2;
L0 = 0.01; %1e-5;

x1 = linspace(0,L1,N1+1); % mesh in x1
x2 = linspace(0,L2,N2+1); % mesh in x2
[xx1m,xx2m] = meshgrid(x1,x2); % all nodes

xx1 = xx1m(:); % x1 coordinates in a vvector
xx2 = xx2m(:); % x2 coordinates in a vvector

Nuzel = (N1+1)*(N2+1) % number of nodes
Nuzelp = N1*N2; % number of nodes without "periodic"

elem = MyDelaunay(N1,N2); % set of elements (triangles) 
Nelem = size(elem,1); % number of elements
% KresliElem(elem,xx1,xx2,Nelem); % drawing

% denoting boundary and "periodic" nodes:
% Nuzel1 = Nuzel;
uzel_hra = Najdi_uzel_hra(xx1,xx2,Nuzel,L1,L2); % find boundary nodes
% Nuzel1 = Nuzel-sum(uzel_hra);

% building matrix A and RHS B for all nodes:
A = zeros(Nuzel);
Apre  = A;
B = zeros(Nuzel,1);
% go along elements:
for k = 1:Nelem
    cu = elem(k,:); % cisla uzlu
    A(cu,cu) = A(cu,cu)+Integruj_a(k,elem,xx1,xx2,N1,N2);
    Apre(cu,cu) = Apre(cu,cu)+Integruj_apre(k,elem,xx1,xx2,N1,N2);
    B(cu) = B(cu)+Integruj_f(k,elem,xx1,xx2,N1,N2,gc_mac);
end;

A3 = A;
A3pre = Apre;
B3 = -B;
% periodicke OP:
A = PeriodOP_A(A3,N1,N2,uzel_hra,Nuzel);
Apre = PeriodOP_A(A3pre,N1,N2,uzel_hra,Nuzel);
B = PeriodOP_B(B3,N1,N2,uzel_hra,Nuzel);

tic;
[x0,st,ppp,toc1] = CGP_3_without(A,B,B*0,1e-8,max_steps);
time1a = toc;
st1a = st;
Apre*ones(Nuzelp,1);
pom = [Apre,ones(Nuzelp,1);ones(1,Nuzelp),0];
pomi = inv(pom);
pomi = pomi(1:Nuzelp,1:Nuzelp);

fpom = Apre(:,1);
fpomf = real(fft2(reshape(fpom,N2,N1)));
fpomf(1,1) = 1;
fpomfi = 1./fpomf;
fpomfi(1,1) = 0;

tic;
[x0,st,ppp,toc1] = CGP_3_pre(A,B,B*0,1e-8,max_steps,pomi);
time1b = toc;
st1b = st;
tic;
[x0,st,ppp,toc1] = CGP_3_pre_fft(sparse(A),B,B*0,1e-8,max_steps,N1,N2,fpomfi);
time1c = toc;
st1c = st;
steps_diffusion = [st1a,st1b,st1c]
time_diffusion = [time1a,time1b,time1c]


% prolong solution to "periodic" nodes:
UU = Res_period(x0,uzel_hra,Nuzel,N1,N2);
U1 = reshape(UU,N2+1,N1+1);

subplot(2,3,1);
cla; hold on;
surf(xx1m,xx2m,U1,'FaceAlpha',0.5);
title('\phi');

gcx = gc_mac(1)*xx1+gc_mac(2)*xx2; % solution including macro-gradient
U2 = UU+gcx;

% constants beta:
betasol = 0;
betapor = 0;
kk = 0;
for k = 1:Nuzel
    x1p = xx1(k);
    x2p = xx2(k);
    apom = Fcea(x1p,x2p,N1,N2);
    if (norm(apom)<0.5) betasol = betasol+U2(k); kk = kk+1;
    else betapor = betapor +U2(k); end;  %-norm(gc)/2;
end;
betasol = -betasol/kk
betapor = -betapor/(Nuzel-kk)
eta = (Nuzel-kk)/Nuzel

for k = 1:Nuzel
    x1p = xx1(k);
    x2p = xx2(k);
    apom = Fcea(x1p,x2p,N1,N2);
    if (norm(apom)<.5) U2(k) = U2(k)+betasol; % solid
    else U2(k) = U2(k)+betapor+c_mac/eta; % pores
    end;
end;
CC = U2;
subplot(2,3,2);
U2 = reshape(U2,N2+1,N1+1); % concentration c(x) formula (12)
surf(xx1m,xx2m,U2,'FaceAlpha',0.5);
title('\phi+\nabla c_{mac}\cdot x+\beta');
view(0,90)

% eps_eig:
eigstrain_elem = zeros(Nelem,1);
eigstrain_uzel = zeros(Nuzel,1);
for k = 1:Nelem
    c = elem(k,:);
    x1p = xx1(c);
    x2p = xx2(c); 
    x1s = mean(x1p);
    x2s = mean(x2p);
    upom = CC(c);
    apom = Fcea(x1s,x2s,N1,N2);
    precipi = max([0,sum(upom)/3-threshold_c]);
    moles_precipi = precipi*h1*h2;
    volume_precipi = moles_precipi*molar_mass/density;    
    eigen_strain = max([0,(volume_precipi-h1*h2)/h1/h2]);
    precipi_fraction = min([0.9,abs(volume_precipi)/h1/h2]);
    if (norm(upom)>0.5) precipi_fraction = 0; end;
    porosity_reduction = 1-precipi_fraction;    
    eigstrain_elem(k) = eigen_strain;
    eigstrain_uzel(c) = eigen_strain;
end;

subplot(2,3,3);
surf(xx1m,xx2m,reshape(eigstrain_uzel,N2+1,N1+1));
title('\epsilon_{eig} ');
view(0,90)

% stress development:
% building matrix K and RHS KB for all nodes:
K1 = zeros(Nuzel);
K2 = zeros(Nuzel);
K12 = zeros(Nuzel);
K1pre = zeros(Nuzel);
K2pre = zeros(Nuzel);
K12pre = zeros(Nuzel);

KB1 = zeros(Nuzel,1);
KB2 = zeros(Nuzel,1);
% jdu po prvcich:
for k = 1:Nelem
    cu = elem(k,:); % numbers of nodes

    pomc = Integruj_C(k,elem,xx1,xx2,N1,N2);
    K1(cu,cu) = K1(cu,cu)+pomc(1:3,1:3);
    K2(cu,cu) = K2(cu,cu)+pomc(4:6,4:6);
    K12(cu,cu) = K12(cu,cu)+pomc(1:3,4:6);

    pomcp = Integruj_C_pre(k,elem,xx1,xx2,N1,N2);
    K1pre(cu,cu) = K1pre(cu,cu)+pomcp(1:3,1:3);
    K2pre(cu,cu) = K2pre(cu,cu)+pomcp(4:6,4:6);
    K12pre(cu,cu) = K12pre(cu,cu)+pomcp(1:3,4:6);

    % B1(cu) = B1(cu)+Integruj_Kf(k,elem,xx1,xx2,N1,N2,gc);
    % bpom = pomc*(eps_eig_elem(k)*[1;1;0]+eps_mac);
    % KB1(cu) = KB1(cu) + bpom;
    epi = eigstrain_elem(k);
    bpom = Integruj_Kf(k,elem,xx1,xx2,N1,N2,eps_mac,epi);
    KB1(cu) = KB1(cu) + bpom(1:3);
    KB2(cu) = KB2(cu) + bpom(4:6);
end;

K1p = PeriodOP_A(K1,N1,N2,uzel_hra,Nuzel);
K2p = PeriodOP_A(K2,N1,N2,uzel_hra,Nuzel);
K12p = PeriodOP_A(K12,N1,N2,uzel_hra,Nuzel);

K1prep = PeriodOP_A(K1pre,N1,N2,uzel_hra,Nuzel);
K2prep = PeriodOP_A(K2pre,N1,N2,uzel_hra,Nuzel);
K12prep = PeriodOP_A(K12pre,N1,N2,uzel_hra,Nuzel);

B1p = PeriodOP_B(KB1,N1,N2,uzel_hra,Nuzel);
B2p = PeriodOP_B(KB2,N1,N2,uzel_hra,Nuzel);
BB = [B1p;B2p];

KK = [K1p,K12p;K12p',K2p];
KKpre = [K1prep,K12prep;K12prep',K2prep];

% U = pinv(KK)*BB;

tic;
[x0,st,ppp,toc1] = CGP_3_without(KK,BB,BB*0,1e-8,max_steps);
time2c = toc;
st2a = st;
nn1 = norm(KK*x0-BB)

fpom1 = K1prep(:,1);
fpom1f = real(fft2(reshape(fpom1,N2,N1)));
fpom1f(1,1) = 1;

fpom2 = K2prep(:,1);
fpom2f = real(fft2(reshape(fpom2,N2,N1)));
fpom2f(1,1) = 1;

fpom12 = K12prep(:,1);
fpom12f = real(fft2(reshape(fpom12,N2,N1)));

fpom21 = K12prep(1,:)';
fpom21f = real(fft2(reshape(fpom21,N2,N1)));

for j2 = 1:N1
    for j1 = 1:N2
        p1 = [fpom1f(j1,j2),fpom12f(j1,j2);fpom21f(j1,j2),fpom2f(j1,j2)];
        p1i = inv(p1);
        fpom1f(j1,j2) = p1i(1,1);
        fpom12f(j1,j2) = p1i(1,2);
        fpom21f(j1,j2) = p1i(2,1);
        fpom2f(j1,j2) = p1i(2,2);
    end;
end;

% ker = [ones(Nuzelp,1),zeros(Nuzelp,1),-xx2p;zeros(Nuzelp,1),ones(Nuzelp,1),xx1p];
ker2 = [ones(Nuzelp,1),zeros(Nuzelp,1);zeros(Nuzelp,1),ones(Nuzelp,1)];

% pom = [KKpre,ker;ker',zeros(3)];
pom2 = [KKpre,ker2;ker2',zeros(2)];

pomi = inv(pom2);
pomi = pomi(1:2*Nuzelp,1:2*Nuzelp);
tic;
[x0,st,ppp,toc1] = CGP_3_pre(KK,BB,BB*0,1e-8,max_steps,pomi);
time2a = toc;
st2b = st;

tic;
[x0,st,ppp,toc1] = CGP_3_pre_elast_fft(sparse(KK),BB,BB*0,1e-8,max_steps,...
    N1,N2,pomi,fpom1f,fpom12f,fpom2f);
time2b = toc;
st2c = st;
% return
nn2 = norm(KK*x0-BB)
steps_elast = [st2a,st2b,st2c]
time_elast = [time2c,time2a,time2b]


% crack initiation and propagation:
U = x0;
U1 = U(1:N1*N2);
U2 = U(N1*N2+1:end);
U1per = Res_period(U1,uzel_hra,Nuzel,N1,N2);
U2per = Res_period(U2,uzel_hra,Nuzel,N1,N2);
[Psiplus,Psiplus_nodes,epsxx,epsyy] = Psi_plus_elem(U1per,U2per,elem,Nelem,xx1,xx2,N1,N2,eps_mac);
for k = 1:Nelem
    c = elem(k,:);
    x1p = xx1(c);
    x2p = xx2(c); 
    x1s = mean(x1p);
    x2s = mean(x2p);
    upom = Psiplus_nodes(c);
    apom = Fcea(x1s,x2s,N1,N2);
    if (norm(apom)>0.5) Psiplus_nodes(c) = 0; epsxx(c) = 0; end;    
end;

subplot(2,3,4);
surf(reshape(Psiplus_nodes,N2+1,N1+1));
title('Positive strain energy \psi^+');
view(0,90)

subplot(2,3,5);
surf(reshape(epsxx,N2+1,N1+1));
title('Strain \epsilon_{xx}');
view(0,90)

% damage computation:

Adam = zeros(Nuzel);
Adam_pre = zeros(Nuzel);
Bdam = zeros(Nuzel,1);
% jdu po prvcich:
for k = 1:Nelem
    cu = elem(k,:); % cisla uzlu
    Adam(cu,cu) = Adam(cu,cu)+Integruj_Adam(k,elem,xx1,xx2,N1,N2,Gc,L0,Psiplus);
    Adam_pre(cu,cu) = Adam_pre(cu,cu)+Integruj_Adam_pre(k,elem,xx1,xx2,N1,N2,Gc,L0,Psiplus);
    Bdam(cu) = Bdam(cu)+Integruj_Bdam(k,elem,xx1,xx2,N1,N2,Psiplus);
end;
% periodicke OP:
Adamp = PeriodOP_A(Adam,N1,N2,uzel_hra,Nuzel);
Adamp_pre = PeriodOP_A(Adam_pre,N1,N2,uzel_hra,Nuzel);
Bdamp = PeriodOP_B(Bdam,N1,N2,uzel_hra,Nuzel);
Dam = Adamp\Bdamp;
Dam(Dam<0) = 0;
Dam = Res_period(Dam,uzel_hra,Nuzel,N1,N2);
for k = 1:Nuzel
    x1p = xx1(k);
    x2p = xx2(k);
    apom = Fcea(x1p,x2p,N1,N2);
    if (norm(apom)>0.5) Dam(k) = 0; % solid
    end;
end;

[x0,st,ppp,toc1] = CGP_3_without(Adamp,Bdamp,Bdamp*0,1e-8,1000);
st3a = st;
norm(Adamp*x0-Bdamp)
ppA = inv(Adamp_pre);
[x0,st,ppp,toc1] = CGP_3_pre(Adamp,Bdamp,Bdamp*0,1e-8,1000,ppA);
st3b = st;
damage_steps = [st3a,st3b]
norm(Adamp*x0-Bdamp);

subplot(2,3,6);
surf(reshape(Dam,N2+1,N1+1));
title('Damage d');
view(0,90);

return

% KresliReseni(elem,xx1,xx2,Nelem,upom); % draw the last solution
% save('MaticeMKP1',"A","B");

%====================================================================
%====================================================================
%====================================================================

function z = Fcea(x1,x2,N1,N2)

% k1 = round(x1*N1+0.5);
% k2 = round(x2*N2+0.5);
% z = oblast(k1,k2)*eye(2);
cc = 1 ;
cc2 = cc^2;
z =1 ; % diffusion coef. pores
if ((x1-0.35*cc)^2+(x2-0.35*cc)^2<0.05*cc2 ||...
        (x1-0.65*cc)^2+(x2-0.65*cc)^2<0.05*cc2)
    z = 1e-6; % diffusion coef. solid
end;
z = z*eye(2);

function z = Fcea_pre(x1,x2,N1,N2)
z = eye(2);

function [E,nu] = FceEaNu(x1,x2,N1,N2)
E = 1;  % solid 10, pores 1; 
nu = 0.45;  % solid 0.2, pores 0.45
cc = 1 ;
cc2 = cc^2;
if ((x1-0.35*cc)^2+(x2-0.35*cc)^2<0.05*cc2 ||...
        (x1-0.65*cc)^2+(x2-0.65*cc)^2<0.05*cc2)
    E = 10;
    nu = 0.2;
end;

function [E,nu] = FceEaNu_pre(x1,x2,N1,N2)
E = 1;  % solid 10, pores 1; 
nu = 0.2;  % solid 0.2, pores 0.45

function z = FcegD(x1,x2)
z = 0;

function z = KresliElem(elem,xx1,xx2,Nelem)
cla; hold on;
for k = 1:Nelem
    t = elem(k,:);
    plot(xx1([t,t(1)]),xx2([t,t(1)]));
end;

function z = Najdi_uzel_hra(xx1,xx2,Nuzel,L1,L2)
epsi = 1e-8;
cc = 0; 
uzel_hra = zeros(Nuzel,1);
for k = 1:Nuzel
    if (abs(xx1(k)-L1)<epsi || abs(xx1(k)-L1*cc)<epsi || ...
            abs(xx2(k)-L2)<epsi || abs(xx2(k))<epsi)
        uzel_hra(k) = 1;        
    end;
    if (abs(xx2(k)-L2)<epsi )
        uzel_hra(k) = 2;        
    end; 
    if (abs(xx1(k)-L1)<epsi)
        uzel_hra(k) = 3;        
    end;   
   
end;
z = uzel_hra;

function z = KresliReseni(elem,xx1,xx2,Nelem,u)
cla; hold on;
for k = 1:Nelem
    t = elem(k,:);
    surf([xx1(t(1)),xx1(t(2));xx1(t(3)),xx1(t(3))],...
        [xx2(t(1)),xx2(t(2));xx2(t(3)),xx2(t(3))],...
        [u(t(1)),u(t(2));u(t(3)),u(t(3))],'FaceAlpha',0.5);
end;
view(10,15);

function z = Integruj_a(k,elem,xx1,xx2,N1,N2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
z = pom'*Fcea(x1c,x2c,N1,N2)*pom*d;

function z = Integruj_apre(k,elem,xx1,xx2,N1,N2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
z = pom'*Fcea_pre(x1c,x2c,N1,N2)*pom*d;

function z = Integruj_Adam(k,elem,xx1,xx2,N1,N2,Gc,L0,Hplus)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
p1 = Gc*L0/2*pom'*pom*d;
p2 = (Gc/2/L0+Hplus(k))*[2,1,1;1,2,1;1,1,2]*d/12;
z = p1+p2;

function z = Integruj_Adam_pre(k,elem,xx1,xx2,N1,N2,Gc,L0,Hplus)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
p1 = Gc*L0/2*pom'*pom*d;
p2 = (Gc/2/L0)*[2,1,1;1,2,1;1,1,2]*d/12;
z = p1+p2;

function z = Integruj_Bdam(k,elem,xx1,xx2,N1,N2,Hplus);
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
z = Hplus(k)*ones(3,1)*d/2;

function z = Integruj_C(k,elem,xx1,xx2,N1,N2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
[fc,nu] = FceEaNu(x1c,x2c,N1,N2);
b = (1-2*nu)/2;
pom;
z1 = fc*pom'*[1-nu,0;0,b]*pom*d;
z2 = fc*pom'*[b,0;0,1-nu]*pom*d;
z12 = fc*pom'*[0,nu;b,0]*pom*d;
z = [z1,z12;z12',z2];

function z = Integruj_C_pre(k,elem,xx1,xx2,N1,N2)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
[fc,nu] = FceEaNu_pre(x1c,x2c,N1,N2);
b = (1-2*nu)/2;
pom;
z1 = fc*pom'*[1-nu,0;0,b]*pom*d;
z2 = fc*pom'*[b,0;0,1-nu]*pom*d;
z12 = fc*pom'*[0,nu;b,0]*pom*d;
z = [z1,z12;z12',z2];

function z = Integruj_f(k,elem,xx1,xx2,N1,N2,gc)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
z = pom'*Fcea(x1c,x2c,N1,N2)*gc*d;

function z = Integruj_Kf(k,elem,xx1,xx2,N1,N2,eps_mac,eps_eig)
t = elem(k,:);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
[fc,nu] = FceEaNu(x1c,x2c,N1,N2);
b = (1-2*nu)/2;
pomv = -eps_mac+eps_eig;
pomv = [pomv(1);0;0;pomv(2)];
z1 = fc*pom'*[1-nu,0;0,b]*pomv(1:2)*d;
z2 = fc*pom'*[b,0;0,1-nu]*pomv(3:4)*d;
z12 = fc*pom'*[0,nu;b,0]*pomv(3:4)*d;
z21 = fc*pom'*[0,b;nu,0]*pomv(1:2)*d;
z = [z1+z12;z21+z2]; 

function z = Res_period(U,uzel_hra,Nuzel,N1,N2)
% prolong solution to "periodic" nodes:
UUper = zeros(Nuzel,1);
UUper(uzel_hra<=1,:) = U;
for r = 1:Nuzel %Nuzel:-1:1
    if (uzel_hra(r)==2) 
        UUper(r) = UUper(r-N2);        
    end;
    if (uzel_hra(r)==3) 
        UUper(r) = UUper(r-N1*(N2+1));        
    end;    
end;
z = UUper;

function A = PeriodOP_A(AA,N1,N2,uzel_hra,Nuzel)
A = AA;
for r = Nuzel:-1:1
    if (uzel_hra(r)==3) 
        A(:,r-N1*(N2+1)) = A(:,r-N1*(N2+1))+A(:,r);
        A(r-N1*(N2+1),:) = A(r-N1*(N2+1),:)+A(r,:);
        A(:,r) = [];
        A(r,:) = [];                
    end;
    if (uzel_hra(r)==2) 
        A(:,r-N2) = A(:,r-N2)+A(:,r);
        A(r-N2,:) = A(r-N2,:)+A(r,:);
        A(:,r) = [];
        A(r,:) = [];              
    end;    
end;
    
function B = PeriodOP_B(BB,N1,N2,uzel_hra,Nuzel)
B = BB;
for r = Nuzel:-1:1
    if (uzel_hra(r)==3)         
        B(r-N1*(N2+1),:) = B(r-N1*(N2+1),:)+B(r,:);
        B(r,:) = [];              
    end;
    if (uzel_hra(r)==2)         
        B(r-N2,:) = B(r-N2,:)+B(r,:);
        B(r,:) = [];        
    end;    
end;

function [zel,znod,epsxx,epsyy] = Psi_plus_elem(U1per,U2per,elem,Nelem,xx1,xx2,N1,N2,eps_mac)
% lambda == nu   a   mi == (1-2*nu)/2
zel = zeros(Nelem,1);
znod = zeros((N1+1)*(N2+1),1);
epsxx = zeros((N1+1)*(N2+1),1);
epsyy = zeros((N1+1)*(N2+1),1);
for k = 1:Nelem
t = elem(k,:);
u1 = U1per(t);
u2 = U2per(t);
x1c = sum(xx1(t)/3);
x2c = sum(xx2(t)/3);
pom  = [xx1(t),xx2(t),[1;1;1]];
d = [pom(1,1:2)-pom(3,1:2);pom(2,1:2)-pom(3,1:2)];
d = abs(det(d))/2;
pom = inv(pom);
pom = pom(1:2,:);
[fc,nu] = FceEaNu(x1c,x2c,N1,N2);

du1 = pom*u1; 
du2 = pom*u2;
du1(1) = du1(1) + eps_mac(1);
du2(2) = du2(2) + eps_mac(2);
epsxx(t) = epsxx(t) +  du1(1);
epsyy(t) = epsyy(t) +  du2(2);

mu1 = fc/2/(1+nu);
lam1 = fc*nu/(1+nu)/(1-2*nu);

tr_eps_2 = 0.25*max(0,(du1(1)+du2(2)+abs(du1(1)+du2(2))))^2; % /4 OK
% tr_eps_2 = 0.25*(du1(1)+du2(2)+abs(du1(1)+du2(2)))^2;
pom = [du1,du2];
eps_dev = pom - (pom(1,1)+pom(2,2))/2*eye(2); % OK
eps_dev = eps_dev(:);
eps_dev = sum(eps_dev.^2);
p = 0.5*(lam1+2*mu1/2)*tr_eps_2+mu1*eps_dev;
zel(k) = p;
znod(t) = znod(t)+zel(k);
end;
znod = znod/6;
epsxx = epsxx/6;
epsyy = epsyy/6;

function el = MyDelaunay(N1,N2)
el = zeros(N1*N2*2,3);
kde = 1;
for k1 = 1:N1
    d = (k1-1)*(N2+1)+1;
    for k2 = 1:N2
        el(2*kde-1,:) = [d,d+1,d+N2+1];
        el(2*kde,:) = [d+1,d+N2+1,d+N2+2];
        kde = kde+1;
        d = d+1;
    end;
end;

function pz = Au_via_fft(fpomfi,z,N1,N2)
% fpomf = real(fft2(reshape(fpom,N2,N1)))
% fpomf(1,1) = 1;
% fpomfi = 1./fpomf;
% fpomfi(1,1) = 0;
zf = fft2(z,N2,N1);
pzf = fpomfi.*zf;
pz = ifft2(pzf,N2,N1);
pz  = pz(:);


