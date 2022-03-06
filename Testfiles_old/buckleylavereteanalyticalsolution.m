%buckley leveratt analytical solution

muw = 1e-3;muo = 3e-3;
sw = 1:-0.0001:0.0;
x = 0:0.0001:1; 
dfractional = @(s) ((2*s./muw).*(s.^2./muw + (1-s).^2./muo) - (s.^2/muw).*(2*s/muw - 2*s.*(1-s)./muo))./(s.^2/muw + (1-s).^2./muo).^2;
%fractiona = 
xd = @(td,sw) td.*dfractional(sw);
td = 1;
i = sw>0.605;
xdd = xd(td,sw(i));
x(i) = xd(td,sw(i));
dx = (1-xdd(end))/(length(sw(~i))-1);
sw(~i) = 0;

x(~i) = xdd(end):dx:1;
% x(ii) = xdd(end);
plot(x,sw,'--','linewidth',2,'color','r')
xlabel('non-dimensionalized distance')
ylabel('saturation,Sw')
legend('Buckley-leverett analytical solution')