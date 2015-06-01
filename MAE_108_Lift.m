clear all;
clc;
format long;

%X Z Position of Clark Y-14 
XClarkTop = xlsread('DataCP.xls','XY Coor','B4:B64');
ZClarkTop = xlsread('DataCP.xls','XY Coor','C4:C64');
XClarkBot = xlsread('DataCP.xls','XY Coor','B66:B126');
ZClarkBot = xlsread('DataCP.xls','XY Coor','C66:C126');

%Data Acquistion Points (through extrapolation)
%Initialize Points
Xtop = [0 .075 .1 .2 .3 .4 .5 .6 .7 .8];
Xbot = [.075 .1 .2 .3 .4 .5 .6 .7];

%Extrapolation through the Clark Coordinates
Ztop = interp1(XClarkTop,ZClarkTop,Xtop,'spline');
Zbot = interp1(XClarkBot,ZClarkBot,Xbot,'spline');

%Top Pressure Coefficient at angles of attack 0, 4, 8, and 20
TCp0 = fliplr(xlsread('DataCP.xls','Pressure','C10:C19'));
TCp4 = fliplr(xlsread('DataCP.xls','Pressure','F10:F19'));
TCp8 = fliplr(xlsread('DataCP.xls','Pressure','I10:I19'));
TCp20 = fliplr(xlsread('DataCP.xls','Pressure','L10:L19'));

%Bottom Pressure Coefficient at angles of attack 0, 4, 8 , and 20
BCp0 = xlsread('DataCP.xls','Pressure','C20:C27');
BCp4 = xlsread('DataCP.xls','Pressure','F20:F27');
BCp8 = xlsread('DataCP.xls','Pressure','I20:I27');
BCp20 = xlsread('DataCP.xls','Pressure','L20:L27');

%Plots
figure(1)
plot(Xtop, Ztop,'x', Xbot, Zbot,'x')
title('NACA 4412 AIRFOIL DATA POSITION')
figure(2)
plot(Xtop, TCp0,'x', Xbot, BCp0,'x')
title('PRESSURE DISTRIBUTION AS A FUNCTION OF CHORD DISTANCE, 0 DEGREES ANGLE OF ATTACK')
figure(3)
plot(Xtop, TCp4,'x', Xbot, BCp4,'x')
title('PRESSURE DISTRIBUTION AS A FUNCTION OF CHORD DISTANCE, 4 DEGREES ANGLE OF ATTACK')
figure(4)
plot(Xtop, TCp8,'x', Xbot, BCp8,'x')
title('PRESSURE DISTRIBUTION AS A FUNCTION OF CHORD DISTANCE, 8 DEGREES ANGLE OF ATTACK')
figure(5)
plot(Xtop, TCp20,'x', Xbot, BCp20,'x')
title('PRESSURE DISTRIBUTION AS A FUNCTION OF CHORD DISTANCE, 20 DEGREES ANGLE OF ATTACK')

%Resultant of pressure in the direction normal to the chord
dx1 = 0*Xtop;
dzdx1 = 0*Xtop;
for i=1:length(Xtop)
    if i == 1
        dx1(i) = (Xtop(2) + Xtop(1))/2;
        dzdx1(i) = (Ztop(2)-Ztop(1))/(Xtop(2)-Xtop(1));
    end;
    if i < length(Xtop) && i > 1
        dx1(i) = (Xtop(i+1)+Xtop(i))/2 - (Xtop(i)+Xtop(i-1))/2;
        dzdx1(i) = 1/((Xtop(i+1)-Xtop(i))*(Xtop(i)-Xtop(i-1))*((Xtop(i)-Xtop(i-1))+(Xtop(i+1)-Xtop(i))))*(Ztop(i+1)*(Xtop(i)-Xtop(i-1))^2-Ztop(i-1)*(Xtop(i+1)-Xtop(i))^2+Ztop(i)*((Xtop(i+1)-Xtop(i))^2-(Xtop(i)-Xtop(i-1))^2));
    end;
    if i == length(Xtop)
        dx1(i) = Xtop(i) - (Xtop(i)+Xtop(i-1))/2;
        dzdx1(i) = (Ztop(i)-Ztop(i-1))/(Xtop(i)-Xtop(i-1));
    end;
end;
dx2 = 0*Xbot;
dzdx2 = 0*Xbot;
for i=1:length(Xbot)
    if i == 1
        dx2(i) = (Xbot(2) + Xbot(1))/2;
        dzdx2(i) = (Zbot(2)-Zbot(1))/(Xbot(2)-Xbot(1));
    end;
    if i < length(Xbot) && i > 1
        dx2(i) = (Xbot(i+1)+Xbot(i))/2 - (Xbot(i)+Xbot(i-1))/2;
        dzdx2(i) = 1/((Xbot(i+1)-Xbot(i))*(Xbot(i)-Xbot(i-1))*((Xbot(i)-Xbot(i-1))+(Xbot(i+1)-Xbot(i))))*(Zbot(i+1)*(Xbot(i)-Xbot(i-1))^2-Zbot(i-1)*(Xbot(i+1)-Xbot(i))^2+Zbot(i)*((Xbot(i+1)-Xbot(i))^2-(Xbot(i)-Xbot(i-1))^2));
    end;
    if i == length(Xbot)
        dx2(i) = Xtop(i) - (Xbot(i)+Xbot(i-1))/2;
        dzdx2(i) = (Zbot(i)-Zbot(i-1))/(Xbot(i)-Xbot(i-1));
    end;
end;

%Check if values are correct
fprintf('Values of dx for the bottom integral = \n');
disp(dx1);
fprintf('\n\nSum of dx for the bottom integral = %d \n\n', sum(dx1));
fprintf('Values of dx for the top integral = \n');
disp(dx2);
fprintf('\nSum of dx for the top integral = %d \n\n', sum(dx2));

%0 Angle of Attack Pressure Normal
Pn0TC=0*Xtop;
for i=1:length(Xtop)
    Pn0TC(i) = dx1(i)*TCp0(i);
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 0 ANGLE OF ATTACK, BOTTOM CONTRIBUTION = %d \n\n', sum(Pn0TC));

Pn0BC=0*Xbot;
for i=1:length(Xbot)
    Pn0BC(i) = -(dx2(i)*BCp0(i));
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 0 ANGLE OF ATTACK, TOP CONTRIBUTION = %d \n\n', sum(Pn0BC));
fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORAML TO THE CHORD, \n 0 ANGLE OF ATTACK: BOTTOM CONTRIBUTION - TOP CONTRIBUTION = %d \n\n\n',sum(Pn0BC)+sum(Pn0TC));

%4 Angle of Attack Pressure Normal
Pn4TC=0*Xtop;
for i=1:length(Xtop)
    Pn4TC(i) = dx1(i)*TCp4(i);
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 4 ANGLE OF ATTACK, BOTTOM CONTRIBUTION = %d \n\n', sum(Pn4TC));

Pn4BC=0*Xbot;
for i=1:length(Xbot)
    Pn4BC(i) = -(dx2(i)*BCp4(i));
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 4 ANGLE OF ATTACK, TOP CONTRIBUTION = %d \n\n', sum(Pn4BC));
fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORAML TO THE CHORD, \n 4 ANGLE OF ATTACK: BOTTOM CONTRIBUTION - TOP CONTRIBUTION = %d \n\n\n',sum(Pn4BC)+sum(Pn4TC));

%8 Angle of Attack Pressure Normal
Pn8TC=0*Xtop;
for i=1:length(Xtop)
    Pn8TC(i) = dx1(i)*TCp8(i);
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 8 ANGLE OF ATTACK, BOTTOM CONTRIBUTION = %d \n\n', sum(Pn8TC));

Pn8BC=0*Xbot;
for i=1:length(Xbot)
    Pn8BC(i) = -(dx2(i)*BCp8(i));
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 8 ANGLE OF ATTACK, TOP CONTRIBUTION = %d \n\n', sum(Pn8BC));
fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORAML TO THE CHORD, \n 8 ANGLE OF ATTACK: BOTTOM CONTRIBUTION - TOP CONTRIBUTION = %d \n\n\n',sum(Pn8BC)+sum(Pn8TC));

%20 Angle of Attack Pressure Normal
Pn20TC=0*Xtop;
for i=1:length(Xtop)
    Pn20TC(i) = dx1(i)*TCp20(i);
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 20 ANGLE OF ATTACK, BOTTOM CONTRIBUTION = %d \n\n', sum(Pn20TC));

Pn20BC=0*Xbot;
for i=1:length(Xbot)
    Pn20BC(i) = -(dx2(i)*BCp20(i));
end;

fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORMAL TO THE CHORD, \n 20 ANGLE OF ATTACK, TOP CONTRIBUTION = %d \n\n', sum(Pn20BC));
fprintf('RESULTANT OF PRESSURE IN THE DIRECTION NORAML TO THE CHORD, \n 20 ANGLE OF ATTACK: BOTTOM CONTRIBUTION - TOP CONTRIBUTION = %d \n\n\n',sum(Pn20BC)+sum(Pn20TC));

%Local Angle calculation
xq = 0:.001:1;
BottomProfileF = interp1(Xtop,Ztop,xq,'spline');
figure(6)
plot(Xtop,Ztop,'o',xq,BottomProfileF)
TopAngle=atan(dzdx1);
% TAClark = atan(ZClarkTop./(XClarkTop));
% TopAngle = interp1(XClarkTop,TAClark,Xtop,'spline');
% for i = 6:10
%     Temp = TopAngle(i);
%     TopAngle(i) = -Temp;
% end;

disp(dzdx2);
TopProfileF = interp1(Xbot,Zbot,xq,'spline');
figure(7)
plot(Xbot,Zbot,'o',xq,TopProfileF)
BottomAngle=atan(dzdx2);
% BAClark = atan(ZClarkBot./(XClarkBot));
% BottomAngle = interp1(XClarkBot,BAClark,Xbot,'spline');
% for i = 3:8
%     Temp = BottomAngle(i);
%     BottomAngle(i) = -Temp;
% end;

%0 Angle of Attack Pressure Parallel
Pp0TC=0*dx1;
for i=1:length(dx1)
    Pp0TC(i)=-dx1(i)*tan(TopAngle(i))*TCp0(i);
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 0 ANGLE OF ATTACK, BOTTOM CONTRIBUTION, = %d \n\n',sum(Pp0TC));

Pp0BC=0*dx2;
for i=1:length(dx2)
    Pp0BC(i)=dx2(i)/cos(BottomAngle(i))*(BCp0(i)*sin(BottomAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 0 ANGLE OF ATTACK, TOP CONTRIBUTION, = %d \n\n',sum(Pp0BC));
fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 0 ANGLE OF ATTACK = %d \n\n\n',sum(Pp0TC)+sum(Pp0BC));

%4 Angle of Attack Pressure Parallel
Pp4TC=0*dx1;
for i=1:length(dx1)
    Pp4TC(i)=-dx1(i)/cos(TopAngle(i))*(TCp4(i)*sin(TopAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 4 ANGLE OF ATTACK, BOTTOM CONTRIBUTION, = %d \n\n',sum(Pp4TC));

Pp4BC=0*dx2;
for i=1:length(dx2)
    Pp4BC(i)=dx2(i)/cos(BottomAngle(i))*(BCp4(i)*sin(BottomAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 4 ANGLE OF ATTACK, TOP CONTRIBUTION, = %d \n\n',sum(Pp4BC));
fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 4 ANGLE OF ATTACK = %d \n\n\n',sum(Pp4BC)+sum(Pp4TC));

%8 Angle of Attack Pressure Parallel
Pp8TC=0*dx1;
for i=1:length(dx1)
    Pp8TC(i)=-dx1(i)/cos(TopAngle(i))*(TCp8(i)*sin(TopAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 8 ANGLE OF ATTACK, BOTTOM CONTRIBUTION, = %d \n\n',sum(Pp8TC));

Pp8BC=0*dx2;
for i=1:length(dx2)
    Pp8BC(i)=dx2(i)/cos(BottomAngle(i))*(BCp8(i)*sin(BottomAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 8 ANGLE OF ATTACK, TOP CONTRIBUTION, = %d \n\n',sum(Pp8BC));
fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 8 ANGLE OF ATTACK = %d \n\n\n',sum(Pp8BC)+sum(Pp8TC));

%20 Angle of Attack Pressure Parallel
Pp20TC=0*dx1;
for i=1:length(dx1)
    Pp20TC(i)=-dx1(i)/cos(TopAngle(i))*(TCp20(i)*sin(TopAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 20 ANGLE OF ATTACK, BOTTOM CONTRIBUTION, = %d \n\n',sum(Pp20TC));

Pp20BC=0*dx2;
for i=1:length(dx2)
    Pp20BC(i)=dx2(i)/cos(BottomAngle(i))*(BCp20(i)*sin(BottomAngle(i)));
end;

fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 20 ANGLE OF ATTACK, TOP CONTRIBUTION, = %d \n\n',sum(Pp20BC));
fprintf('RESULTANT OF PRESSURE PARALLEL TO THE CHORD, \n 20 ANGLE OF ATTACK = %d \n\n\n',sum(Pp20BC)+sum(Pp20TC));

% Drag/Lift 0 Angle of Attack
DragAA0 = (sum(Pn0TC)+sum(Pn0BC))*sin(0*2*pi/360)+(sum(Pp0TC)+sum(Pp0BC))*cos(0*2*pi/360);
LiftAA0 = (sum(Pn0TC)+sum(Pn0BC))*cos(0*2*pi/360)-(sum(Pp0TC)+sum(Pp0BC))*sin(0*2*pi/360);

% Drag/Lift 4 Angle of Attack
DragAA4 = (sum(Pn4TC)+sum(Pn4BC))*sin(4*2*pi/360)+(sum(Pp4TC)+sum(Pp4BC))*cos(4*2*pi/360);
LiftAA4 = (sum(Pn4TC)+sum(Pn4BC))*cos(4*2*pi/360)-(sum(Pp4TC)+sum(Pp4BC))*sin(4*2*pi/360);

% Drag/Lift 8 Angle of Attack
DragAA8 = (sum(Pn8TC)+sum(Pn8BC))*sin(8*2*pi/360)+(sum(Pp8TC)+sum(Pp8BC))*cos(8*2*pi/360);
LiftAA8 = (sum(Pn8TC)+sum(Pn8BC))*cos(8*2*pi/360)-(sum(Pp8TC)+sum(Pp8BC))*sin(8*2*pi/360);

% Drag/Lift 20 Angle of Attack
DragAA20 = (sum(Pn20TC)+sum(Pn20BC))*sin(-4*2*pi/360)+(sum(Pp20TC)+sum(Pp20BC))*cos(20*2*pi/360);
LiftAA20 = (sum(Pn20TC)+sum(Pn20BC))*cos(-4*2*pi/360)-(sum(Pp20TC)+sum(Pp20BC))*sin(20*2*pi/360);

DragCoef = [-DragAA0 -DragAA4 -DragAA8 DragAA20];
LiftCoef = -1*[LiftAA0 LiftAA4 LiftAA8 LiftAA20];
fprintf('\nDrag Coefficient = ');
disp(DragCoef);
fprintf('Lift Coefficient = ');
disp(LiftCoef);

SectionDrag = DragCoef/0.257;
SectionLift = LiftCoef/0.257;
fprintf('\nSection Drag Coefficient = ');
disp(SectionDrag);
fprintf('Section Lift Coefficient = ');
disp(SectionLift);

%Export Values to a SpreadSheet for improved graphs
xlswrite('Output Data.xlsx',XClarkTop,'Pressure Taps','A3');
xlswrite('Output Data.xlsx',ZClarkTop,'Pressure Taps','B3');
xlswrite('Output Data.xlsx',XClarkBot,'Pressure Taps','C3');
xlswrite('Output Data.xlsx',ZClarkBot,'Pressure Taps','D3');
xlswrite('Output Data.xlsx',Xtop(:),'Pressure Taps','E3');
xlswrite('Output Data.xlsx',Ztop(:),'Pressure Taps','F3');
xlswrite('Output Data.xlsx',Xbot(:),'Pressure Taps','G3');
xlswrite('Output Data.xlsx',Zbot(:),'Pressure Taps','H3');
xlswrite('Output Data.xlsx',Pn0TC(:),'Pressure Taps','I3');
xlswrite('Output Data.xlsx',Pp0TC(:),'Pressure Taps','J3');
xlswrite('Output Data.xlsx',Pn0BC(:),'Pressure Taps','K3');
xlswrite('Output Data.xlsx',Pp0BC(:),'Pressure Taps','L3');
xlswrite('Output Data.xlsx',Pn4TC(:),'Pressure Taps','M3');
xlswrite('Output Data.xlsx',Pp4TC(:),'Pressure Taps','N3');
xlswrite('Output Data.xlsx',Pn4BC(:),'Pressure Taps','O3');
xlswrite('Output Data.xlsx',Pp4BC(:),'Pressure Taps','P3');
xlswrite('Output Data.xlsx',Pn8TC(:),'Pressure Taps','Q3');
xlswrite('Output Data.xlsx',Pp8TC(:),'Pressure Taps','R3');
xlswrite('Output Data.xlsx',Pn8BC(:),'Pressure Taps','S3');
xlswrite('Output Data.xlsx',Pp8BC(:),'Pressure Taps','T3');
xlswrite('Output Data.xlsx',Pn20TC(:),'Pressure Taps','U3');
xlswrite('Output Data.xlsx',Pp20TC(:),'Pressure Taps','V3');
xlswrite('Output Data.xlsx',Pn20BC(:),'Pressure Taps','W3');
xlswrite('Output Data.xlsx',Pp20BC(:),'Pressure Taps','X3');
xlswrite('Output Data.xlsx',SectionDrag(:),'Pressure Taps','H20');
xlswrite('Output Data.xlsx',SectionLift(:),'Pressure Taps','L20');
xlswrite('Output Data.xlsx',DragCoef(:),'Pressure Taps','I20');
xlswrite('Output Data.xlsx',LiftCoef(:),'Pressure Taps','M20');



        