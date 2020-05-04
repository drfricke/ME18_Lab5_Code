% Lab 5

aveDia = 0.00955; %m
kSlope = (137-110)/(400-300); %HT Textbook
k = 119.09; %at 60.5C for Brass
Ac = (pi/4)*(aveDia)^2; %m^02
L = 0.394; %m OA length
P = pi*aveDia; %perimeter (m)
Tb = 60.5; %C at base IR
Tt = 24.0; %C at tip IR
aveSpa = 0.05018; %m
h = 12; % W/(m^2 * K)Free Convection - textbook estimate
TAir = 22; % C
thetaB = (Tb - TAir);

m = sqrt(h*P/(k*Ac));

x = 0:L/8:L; %close to aveSpa


%constant Conductivity of Material (k) for reference
Ts_Convection_constant_k = ((cosh(m*(L-(x)))+(h/(m*k))*sinh(m*(L-(x))))/(cosh(m*L)+(h/(m*k))*sinh(m*L)))*(thetaB) + TAir;
Ts_Adiabatic_constant_k = cosh(m*(L-x))/(cosh(m*L))*thetaB + TAir;
Ts_Prescribed_Temperature_constant_k = (((Tt-TAir)/thetaB)*sinh(m*x)+sinh(m*(L-x)))/(sinh(m*L))*thetaB+TAir;
Ts_Infinite_constant_k = exp(-m*x)*(thetaB) + TAir;


%below attempts to get the best possible results with a changing Conductivity of Material (k) with
%respect to location. This also changes m since m=f(k)


%make some empty arrays
Ts_Convection= zeros(1,length(x));
Ts_Adiabatic= zeros(1,length(x));
Ts_Prescribed_Temperature= zeros(1,length(x));
Ts_Infinite= zeros(1,length(x));


%just the EQ %model_convection = (cosh(m(L-x))+(h/(m*k))*sinh(m*(L-x)))/(cosh(m*L)+(h/(m*k))*sinh(m*L));
for i2=1:length(x)
   Ts_Convection(i2) = ((cosh(m*(L-(x(i2))))+(h/(m*k))*sinh(m*(L-(x(i2)))))/(cosh(m*L)+(h/(m*k))*sinh(m*L)))*(thetaB) + TAir;
   k = kSlope*((Ts_Convection(i2)+273.15)-300)+110;
   m=sqrt(h*P/(k*Ac));
end

disp('Ts_Convection')
disp(Ts_Convection)

%%%%
%reset m to normal values
k = 119.09;
m = sqrt(h*P/(k*Ac)); 


for i3=1:length(x)
   Ts_Adiabatic(i3) = cosh(m*(L-x(i3)))/(cosh(m*L))*thetaB + TAir;
   k = kSlope*((Ts_Adiabatic(i3)+273.15)-300)+110;
   m=sqrt(h*P/(k*Ac));
end

disp('Ts_Adiabatic')
disp(Ts_Adiabatic)

%%%%
%reset m to normal values
k = 119.09;
m = sqrt(h*P/(k*Ac)); 


for i4=1:length(x)
   Ts_Prescribed_Temperature(i4) = (((Tt-TAir)/thetaB)*sinh(m*x(i4))+sinh(m*(L-x(i4))))/(sinh(m*L))*thetaB+TAir;
   k = kSlope*((Ts_Convection(i4)+273.15)-300)+110;
   m=sqrt(h*P/(k*Ac));
end

disp('Ts_Prescribed_Temperature')
disp(Ts_Prescribed_Temperature)

%%%%
%reset m to normal values
k = 119.09;
m = sqrt(h*P/(k*Ac)); 

for i=1:length(x)
   Ts_Infinite(i) = exp(-m*x(i))*(thetaB) + TAir;
   k = kSlope*((Ts_Infinite(i)+273.15)-300)+110;
   m=sqrt(h*P/(k*Ac));
end

disp('Ts_Infinite')
disp(Ts_Infinite)

Data = [60.5 58.2 46.7 39.3 34.5 30.3 28.4 27.7 27.1]

figure(1)
plot(x,Data,'x')
title('Actual Data vs. Length')
xlabel('Average Spacing (mm)')
ylabel('Temperature (C)')

figure(2)
subplot(2,2,1)
plot(x,Data,'x',x,Ts_Adiabatic)
title('Adiabatic Modeling vs. Length')
xlabel('Average Spacing (mm)')
ylabel('Temperature (C)')
legend('Data','Model')

subplot(2,2,2)
plot(x,Data,'x',x,Ts_Convection)
title('Convection Modeling vs. Length')
xlabel('Average Spacing (mm)')
ylabel('Temperature (C)')

subplot(2,2,3)
plot(x,Data,'x',x,Ts_Infinite)
title('Infinite Modeling vs. Length')
xlabel('Average Spacing (mm)')
ylabel('Temperature (C)')

subplot(2,2,4)
plot(x,Data,'x',x,Ts_Prescribed_Temperature)
title('Prescribed Temperature vs. Length')
xlabel('Average Spacing (mm)')
ylabel('Temperature (C)')


figure(1)