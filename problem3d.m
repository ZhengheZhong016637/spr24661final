function problem3d
clear all; close all;
h = 0.05;
k = 0.025;
N = 81;
Nt = 400;
u = uinitial;

for i = 1:Nt
    u = GD(u);    
    if i == 40
        u1 = u;        
    elseif i == 80
        u2 = u;        
    elseif i == 120
        u3 = u;        
    elseif i == 160
        u4 = u;           
    elseif i == 200
        u5 = u;      
    elseif i == 240
        u6 = u;      
    elseif i == 280
        u7 = u;      
    elseif i == 320
        u8 = u;      
    elseif i == 360
        u9 = u;           
    end
end
u10 = u;

figure; hold on
plot(u1,'-o');
plot(u2,'LineWidth',2);
plot(u3,'-o');
plot(u4,'LineWidth',2);
plot(u5,'-o');
plot(u6,'LineWidth',2);
plot(u7,'-o');
plot(u8,'LineWidth',2);
plot(u9,'-o');
plot(u10,'LineWidth',2);
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7','t=8','t=9','t=10')

    function unew = GD(u)
        unew = zeros(N,1);
        unew(1)= u(1)-(k/h)*(F(u(1),u(2))-F(u(N),u(1)));        
        for j = 2:N-1
            unew(j) = u(j)-(k/h)*(F(u(j),u(j+1))-F(u(j-1),u(j)));            
        end
        unew(N) = u(N)-(k/h)*(F(u(N),u(1))-F(u(N-1),u(N)));        
    end       

    function ui = uinitial
        ui = zeros(N,1);
        ui(1:40) = 0.1;
        for j = 41:61
            ui(j) = 0.1+0.8*(j-41)*0.05;
        end
        ui(62:N) = 0.9;
    end

    function nflux = F(ul,ur)        
        if ul<=ur %nflux = minf(u)
           if ul<= -1/sqrt(3) && -1/sqrt(3)<=ur
               nflux = min(min(f(ur),f(ul)),-1/sqrt(3)+1/(sqrt(3)^3));
           else
               nflux = min(f(ul),f(ur));
           end
        else %nflux = maxf(u)
            if ur<= 1/sqrt(3) && 1/sqrt(3)<=ul
                nflux = max(max(f(ur),f(ul)),1/sqrt(3)-1/(sqrt(3)^3));
            else
                nflux = max(f(ur),f(ul));
            end
        end        
    end

    function y = f(x)
        y = x-x^3;
    end
end