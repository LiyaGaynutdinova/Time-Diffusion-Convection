function [x0,st] = CGP_3(A,B,xs,toler,steps)
% conjugate gradients:
% given s.p.d. A, r.h.s. B, toler, steps

r0 = B - A*xs; 
r07 = B - A*xs;
nr0 = norm(r0);
x0 = xs;
if (norm(r07)<toler) st = 0; return; end;
p0 = r0;
for k5 = 1:steps
pom1 = A*p0; 
alfa0 = (r0'*r0)/(p0'*pom1);
x1 = x0 + alfa0*p0;
r1 = r0-alfa0*pom1;
if (norm(r1)/nr0<toler) x0 = x1; st = k5; return; end;
beta0 = (r1'*r1)/(r0'*r0);
p1 = r1 + beta0*p0;
p0 = p1;
r0 = r1;
x0 = x1; 
end;
st = -1;















