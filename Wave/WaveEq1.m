
function Wave
% wave propagation in non-homogeneous media, 1D example

global A;
global alfa;
global L;
global E_konst;
global Rho_konst;
A = 0.001;
alfa = 4;
E_konst = 70.3e9;
Rho_konst = 2700;
clc
all_sols = 0; % all solutions during time interval (yes 1, no 0)
method = 3; %  1-only GEM;  2-CG no precond.;  3-CG with FFT-precond.
L = 2; % length of x-interval
T = 3.919e-4; % length t-interval; recom. 3.919e-4 
T = 4e-4;
Tstop = 0.95e-4; % time end
Tstop = 1.6e-4; % time end
Nx = 500; % 6561 steps in x (recom.)
c0 = 5102.6; % recom. %sqrt(E/rho);
c0 = sqrt(E_konst/Rho_konst);
CFL = L/Nx/c0 % recom.
omega = 5*pi*c0/L;  % L=2
toler = 1e-9; % for CG
max_steps = 1000; % for CG

dt = 8*CFL* 0.5; % recom.
dt = dt*0.5 % different time steps 

Nt  = round(T/dt) % number of time steps
x = linspace(0,L,Nx+1);
dx = L/Nx;
dt = T/Nt;
beta = 0.25; % Newmark constants
gama = 0.5;
% nx0 = round(Nx/2)+1; % index of x0
nx0 = 1; % index for Green function
xx0 = zeros(Nx,1);
xx0(nx0) = 1;
Ft = zeros(Nt+1,1);

K = zeros(Nx+1); % stiffness matrix
Kp = zeros(Nx+1); % stiffness matrix - precond.
Drho = zeros(Nx+1); % mass matrix
Drhop = zeros(Nx+1); % mass matrix - precond.
alfa1 = 1e12; % spectral bounds
alfa2 = 0; % spectral bounds
tic;
% bulding matrices:
for j = 1:Nx
    pom_rho = FceRho(j*dx-dx/2);    
    pom_rhop = FceRhop(j*dx-dx/2);    
    pomE = FceE(j*dx-dx/2);    
    pomEp = FceEp(j*dx-dx/2); 
    K(j:j+1,j:j+1) = K(j:j+1,j:j+1) + [-1,1;1,-1]*pomE/dx;
    Kp(j:j+1,j:j+1) = Kp(j:j+1,j:j+1) + [-1,1;1,-1]*pomEp/dx;
    Drho(j:j+1,j:j+1) = Drho(j:j+1,j:j+1) + [1/3,1/6;1/6,1/3]*dx*pom_rho;
    Drhop(j:j+1,j:j+1) = Drhop(j:j+1,j:j+1) + [1/3,1/6;1/6,1/3]*dx*pom_rhop;
    m = beta*dt^2*pomE*[1,-1;-1,1]/dx+[1/3,1/6;1/6,1/3]*dx*pom_rho;
    mp = beta*dt^2*pomEp*[1,-1;-1,1]/dx+[1/3,1/6;1/6,1/3]*dx*pom_rhop;
    ee = eig(inv(mp)*m);
    if (min(ee)<alfa1) alfa1 = min(ee); end;
    if (max(ee)>alfa2) alfa2 = max(ee); end;      
end;

% periodicity:
K(1,:) = K(1,:) + K(Nx+1,:);
K(:,1) = K(:,1) + K(:,Nx+1);
Kp(1,:) = Kp(1,:) + Kp(Nx+1,:);
Kp(:,1) = Kp(:,1) + Kp(:,Nx+1);
Drho(1,:) = Drho(1,:) + Drho(Nx+1,:);
Drho(:,1) = Drho(:,1) + Drho(:,Nx+1);
Drhop(1,:) = Drhop(1,:) + Drhop(Nx+1,:);
Drhop(:,1) = Drhop(:,1) + Drhop(:,Nx+1);
K = K(1:Nx,1:Nx);
Kp = Kp(1:Nx,1:Nx);
Drho = Drho(1:Nx,1:Nx);
Drhop = Drhop(1:Nx,1:Nx);
Krho = -beta*dt^2*K + Drho;
Mrho = -beta*dt^2*Kp + Drhop;
cas1 = toc;
Kpfi = 1./real(fft(Mrho(:,1))); % inversion Kp fft
tic;
U0 = zeros(Nx,1);
U0d = zeros(Nx,1); % u'
U0dd = zeros(Nx,1); % u''
U0(1:Nx) = FceU0(x(1:Nx));
U0d(1:Nx) = FceV0(x(1:Nx));
% U0dd = K*U0; 
gxx0 = Krho\xx0;
gx0x0 = gxx0(nx0);
cas2 = toc;
tic;
if (all_sols==1) 
    Ukres = zeros(Nt+1,Nx);
    Ukres(1,:) = U0;
end;

tn = 0; % time 
aver_st = 0; % averaged number of steps
for n = 1:Nt
    n;
    tn = tn+dt;    
    bn0 = Drho*(U0+dt*U0d+dt^2*(0.5-beta)*U0dd);   
    if (method == 1)
        Ub = Krho\bn0;
    else if (method == 2)
        [Ub,st,ppp,toc1]  = CGP_3_without(Krho,bn0,bn0*0,toler,max_steps);
        stw = st;
    else 
        [Ub,st,ppp,toc1]  = CGP_3_fft(Krho,bn0,bn0*0,toler,max_steps,Kpfi);
        stp = st;
    end; 
    end;    
    aver_st = aver_st + st;
    Fn = (FceUt(tn)-Ub(nx0))/(beta*dt^2*gx0x0);
    Un = beta*dt^2*Fn*gxx0 + Ub;
    Undd = (1/(beta*dt^2))*(Un-U0-dt*U0d-dt^2*(0.5-beta)*U0dd);
    Und = U0d + dt*(1-gama)*U0dd+dt*gama*Undd;
    if (all_sols==1) Ukres(n+1,:) = Un; end;
    U0 = Un;
    U0d = Und;
    U0dd = Undd;
    Ft(n) = Fn;
    if (tn>Tstop) break; end;       
end;
aver_st = aver_st/(n-1)

% return
cas3 = toc;
time123 = [cas1,cas2,cas3]

% subplot(2,1,2)
cla; hold on;
if (all_sols==1) surf(Ukres); 
else 
    plot(x(1:end-1),U0,'b','LineWidth',2);
    % axis([0,2,-0.1e-3,1.1e-3]);
end;
xlabel('x')
ylabel('u(x)')

% drawing of the exact solution for constant data - only one part of the wave:
c0 = sqrt(E_konst/Rho_konst); % speed of the wave
omega = 5*pi*c0/L;  
P = zeros(Nx,1);
dt1 = dx/c0;  % time step (one interval move)
for t1 = 0:dt1:Tstop
    P(2:Nx) = P(1:Nx-1);
    if (t1<pi/omega) P(1) = A*(t1*(pi/omega-t1))^alfa/((pi/omega)^2/4)^alfa; end;
end;
P(end-1:-1:end-round(Nx/2)) = P(1:round(Nx/2));
plot(x(1:Nx)+dx,P,'r--','LineWidth',2);
axis([0,L,-1.5e-4,10.5e-4])

error_final_L2 = norm(P(1:Nx/2)-U0(3:Nx/2+2))*L/Nx*2 % error

% this is demanding for large Nx:
pom = sort(eig(Krho)); % condtioning of the original problem
kapaKM = pom(end)/pom(1)
pomp = sort(eig(inv(Mrho)*Krho)); % conditioning of the precond. problem
kapa_pre = pomp(end)/pomp(1)
alfa1;
alfa2;
kapa_upper_bound = alfa2/alfa1 % estim. of conditioning of the precond. problem
return


%============================================================
%============================================================
%============================================================
%============================================================


function z = FceRho(x)
% z = 1; % T1
z = 2700;
if (x>0.6 && x<1.4) z = 7850; end;

function z = FceRhop(x)
z = 1; % T1
z = 2700; % T2
% z = 10;

function z = FceE(x)
z = 1; % T1
z = 70.3e9; %
if (x>0.6 && x<1.4) z = 211.4e9; end;

function z = FceEp(x)
z = 1; % T1
z = 70.3e9; % T2

function z = FceU0(x) % starting u, periodic
% z = -sin(2*pi*x);
% z = exp(-abs((x-0.5).*(x-1.5))*100); % T1
z = 0; % T2

function z = FceV0(x) % deriative of u, periodic
z = 0; 

function z = FceUt(t) % evolving of the condition in x0
global A;
global alfa;
global L;
global E_konst;
global Rho_konst;
z = 0; % T1
% A = 0.001;
% alfa = 4;
% L = 2;
c0 = sqrt(E_konst/Rho_konst);
% c0 = 5102.6; %sqrt(E/rho);
omega = 5*pi*c0/L;  % L=2
z = 0;
if (t<pi/omega) z = A*(t*(pi/omega-t))^alfa/((pi/omega)^2/4)^alfa; end;

