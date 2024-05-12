function problem1
close all; clear all;
fsz = 20;
mu = 10^6;
Tmaxa = 200;
Tmaxb = 2*10^6;
h = 10^(-3);
Nstep = ceil(Tmaxa/h);
ua = zeros(2,Nstep+1);
ta = h*(1:(Nstep+1))';
ub = zeros(2,1e5);
ua(:,1) = [2;0];
ub(:,1) = [2;0];
tb = zeros(1e5,1);
tol = 1.0e-14;
atol = 10^(-5);
rtol = 10^(-5);
itermax = 20;

% code for part a
% tic
for j = 1:Nstep
    [ua(:,j+1),~] = DIRK2step(ua(:,j),h,tol,itermax);
end
% T = toc;
% disp(T)

% code for part b
T=0;
j=2;
tic
while T<Tmaxb
    [unew,uhat] = DIRK2step(ub(:,j-1),h,tol,itermax);
    e = h*(unew - uhat);
    rhs = atol+rtol*norm(unew);
    if norm(e) < rhs
        ub(:,j)=unew;
        T = T+h;        
        tb(j)=T;
        h = h*1.3;
        j = j+1;
    else
        h = h/1.3;
    end    
end
cpuT = toc;
disp(cpuT)
ub = ub(:,1:j-1);
disp(j)
tb = tb(1:j-1);

figure;hold on
a1=plot(ta,ua(1,:),'Linewidth',2); M1 = "x, part a";
a2=plot(ta,ua(2,:),'Linewidth',2); M2 = "y, part a";
a3=plot(tb,ub(1,:),'Linewidth',2); M3 = "x, part b";
a4=plot(tb,ub(2,:),'Linewidth',2); M4 = "y, part b";
set(gca,'FontSize',fsz);
xlabel('time','FontSize',fsz);
ylabel('x and y','FontSize',fsz);
legend([a1,a2,a3,a4],[M1,M2,M3,M4]);

% figure;
% plot(u(1,:),u(2,:),'LineWidth',2);
% set(gca,'FontSize',fsz);
% xlabel('x','FontSize',fsz);
% ylabel('y','FontSize',fsz);      
    
    function [unp1,uhat] = DIRK2step(u,h,tol,itermax)        
        gamma = 1.0 - 1.0/sqrt(2);
        k1 = func(u);
        for i = 1 : itermax
            k1 = NewtonIterDIRK2(u,h,k1,gamma);
            aux = u+h*gamma*k1;
            if norm(k1 - func(aux)) < tol
                break
            end
        end
        k2 = k1;
        uaux = u + h*(1-gamma)*k1;
        for i =1 : itermax
            k2 = NewtonIterDIRK2(uaux,h,k2,gamma);
            aux = uaux + h*gamma*k2;
            if norm(k2 - func(aux)) < tol
                break
            end
        end
        unp1 = aux;
        uhat = u+h*k2;
    end

    function dy = func(u)
        dy = zeros(2,1);
        dy(1) = u(2);
        dy(2) = mu*(1-u(1)^2)*u(2)-u(1);
    end
    
    function J = Jac(y)
        J = zeros(2,2);
        J(1,1) = 0;
        J(1,2) = 1;
        J(2,1) = -2*mu*y(1)*y(2)-1;
        J(2,2) = mu*(1-y(1)^2);
    end

    function knew = NewtonIterDIRK2(y,h,k,gamma)    
        aux = y + h*gamma*k;
        F = k - func(aux);
        DF = eye(2) - h*gamma*Jac(aux);
        knew = k - DF\F;
    end
end