function problem2
clear all; close all;
tol = 10^-5;
beta = 2;
nc = 45; %number of fixed points on the circle
phi = (1:nc)'*(2*pi/nc);
cfix1 = [0.2*cos(phi)-0.5,0.2*sin(phi)-0.1];
cfix2 = [0.2*cos(phi)+0.5,0.2*sin(phi)+0.1];
fd = @(p) max(max(drectangle(p,-1,1,-1,1),-dcircle(p,-0.5,-0.1,0.2)),-dcircle(p,0.5,0.1,0.2));
[pts,tri] = distmesh2d(fd,@huniform,0.075,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1;cfix1;cfix2]);

D0 = [];
D1 = [];
for i = 1:size(pts,1)
    if dcircle(pts(i,:),-0.5,-0.1,0.2)<tol
        D0 = [D0 i];
    elseif dcircle(pts(i,:),0.5,0.1,0.2)<tol
        D1 = [D1 i];
    end
end
disp(D1);
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,union(D0,D1));
A = sparse(Npts,Npts);
b = sparse(Npts,1);

for j =1:Ntri            
    A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:))+ stima3(pts(tri(j,:),:));        
end

u=sparse(Npts,1);
u(D1)=1;
b=b-A*u;
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

figure;
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
hold on
axis ij
view(2)

function M = stima3(vertices)
    d = size(vertices,2);
    G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
    M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
    pc = (1/3)*sum(vertices);
    M = B1(pc)*M;
    x1 = vertices(1,1);
    y1 = vertices(1,2);
    x2 = vertices(2,1);
    y2 = vertices(2,2);
    x3 = vertices(3,1);
    y3 = vertices(3,2);
    Atri = 0.5*abs(det([x2-x1,x3-x1;y2-y1,y3-y1]));
    M2 = zeros(3,3);
    M2(1,1) = 2*(y2-y3)*(x3-x2);
    M2(1,2) = (y2-y3)*(x1-x3)+(x3-x2)*(y3-y1);
    M2(1,3) = (y2-y3)*(x2-x1)+(x3-x2)*(y1-y2);
    M2(2,1) = M2(1,2);
    M2(2,2) = 2*(y3-y1)*(x1-x3);
    M2(2,3) = (y3-y1)*(x2-x1)+(x1-x3)*(y1-y2);
    M2(3,1) = M2(1,3);
    M2(3,2) = M2(2,3);
    M2(3,3) = 2*(y1-y2)*(x2-x1);
    M2 = B2(pc)*(Atri/(2*Atri)^2)*M2;
    M = M+M2;
end

    function y = B1(p)
        V = cos(2*pi*p(1))+(p(2)-0.1*sin(pi*p(1)))^2;
        M1 = 1+0.5*cos(pi*p(1));
        y = exp(-beta*V)*M1;
    end

    function y = B2(p)
        V = cos(2*pi*p(1))+(p(2)-0.1*sin(pi*p(1)))^2;
        M2 = 0.5*sin(pi*p(1));
        y = exp(-beta*V)*M2;
    end
end