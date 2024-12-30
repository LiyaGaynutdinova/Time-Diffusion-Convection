function [x0,st,ppp,toc1] = CGP_3(A,B,xs,toler,steps,M)
% conjugate gradients:
% given s.p.d. A, r.h.s. B, toler, steps
toc1 = -1;
r0 = B - A*xs; 
r07 = r0;
z0 = M*(B - A*xs);
nr0 = norm(r0);
ppp(1) = nr0;
x0 = xs;
st = 0;
if (norm(r07)<toler) st = 0; return; end;
p0 = z0;
tic;
for k5 = 1:steps
pom1 = A*p0; 
alfa0 = (r0'*z0)/(p0'*pom1);
x1 = x0 + alfa0*p0;
r1 = r0-alfa0*pom1;
ppp(k5+1) = norm(B-A*x1);
if (norm(r1)/nr0<toler) toc1 = toc; x0 = x1; st = k5; return; end;
z1 = M*r1;  
beta0 = (r1'*z1)/(r0'*z0);
p1 = z1 + beta0*p0;
p0 = p1;
r0 = r1;
x0 = x1;
z0 = z1;
end;
toc1 = toc;
st = -1;















