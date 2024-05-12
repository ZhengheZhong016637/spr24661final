function problem3a
X0 = -2:0.1:2;
figure; hold on
for i=1:20
    x2 = 2*0.97+X0(i);
    plot([X0(i) x2],[0 2],'LineWidth',2)
end
for i = 21:31
    x2 = 2*(1-3*(0.1+0.8*X0(i))^2)+X0(i);
    plot([X0(i) x2],[0 2],'LineWidth',2)
end
for i = 31:41
    x2 = -1.43*2+X0(i);
    plot([X0(i) x2],[0 2],'LineWidth',2)
end
end