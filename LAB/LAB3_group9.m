% Define the altitude range
geopotential_altitude = [0:1:25];
r=6356.76;
% Define temperature profile
temperature = zeros(size(geopotential_altitude)); % Initialize temperature array

% Temperature constants
T_constant = 217.569;
for i=1:26
altitude(i)=(r.*geopotential_altitude(i)./(r-geopotential_altitude(i)));
end
altitude
% Calculate temperatures
for i=1:12
temperature(i) = 300 - 6.5 * geopotential_altitude(i); 
end
for i=13:26
temperature(i) = temperature(12);
end

 Plot the graph
figure;
plot(temperature, altitude);
xlabel('Temperature (K)');
ylabel('Altitude (km)');
title('Temperature vs. Altitude');





R = 287.05; 
g0 = 9.80665; 
T0 = 300; 
P0 = 101325; 
L = 6.5/1000; 
r=6356.76;
for i=1:26
    g(i)=g0./(1+altitude(i)./r).^2
end
for i=1:26
pressure(i) = P0 .*(1 + L .* geopotential_altitude(i).*1000 ./ T0) .^ (-g(i) / (R * L));
end
% Calculate density profile using ideal gas law
density = pressure ./ (R .* temperature);

% Plot the graph
figure;
plot(density,altitude);
xlabel('Density (kg/m^3)');
ylabel(' Altitude (km)');
title('Density vs. Altitude');
grid on;
hold off;

% Task 2
S=11*1.2
ms=220;
me=2*3.5 + 6*1.5 + 14;
mbO=50;
wo=ms+me+mbO;
wo=wo.*ones(1,7);
eff =[8:1:13];
for i=1:6
   p=wo(i+1);
Eg=(0.5*density(1)*20^3*S*(0.03 + 0.7786*(0.75.*wo(i+1)./(0.5*density(1)*20^2*S))^2))*7.5;
Ecl=wo(i+1).*30*sin(30)*333.33;
Ec1(i)=wo(i+1).*56*2000./eff(i);
Ec2(i)=(wo(i+1)-50).*56*1625./eff(i);
Ereq(i)= Eg + Ecl + Ec1(i) + Ec2(i) ;
wb(i)=Ereq(i)./(250*3600);
wo(i+1)=wo(i+1)+wb(i);
j=i;
while wo(j+1)-p < 1e-6
    Eg(i)=(0.5*density(1)*20^3*S*(0.03 + 0.7786*(0.75.*wo(j+1)./(0.5*density(1)*20^2*S))^2))*7.5;
Ecl(i)=wo(j+1).*30*sin(30)*333.33;
Ec1(i)=wo(j+1).*56*2000./eff(i);
Ec2(i)=(wo(j+1)-50).*56*1625./eff(i);
Ereq(i)= Eg + Ecl + Ec1(j) + Ec2(j) ;
wb(i)=Ereq(i)./(250*3600);
p=wo(j+1);
wo(j+1)=wo(j+1)+wb(i);
wo(i+1)=wo(j+1);
end 
end

plot(eff,wb, 'b-', 'LineWidth', 2);
xlabel('L/D ');
ylabel(' Weight of the battery');
title('L/D vs. wb');
 
plot(eff,wo(2:7), 'b-', 'LineWidth', 2);
xlabel('L/D');
ylabel(' Total weight');
title('L/D vs. wo');




    













