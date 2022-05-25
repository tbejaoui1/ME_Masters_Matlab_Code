
Position= readtable('Position_mobile.xlsx');

lla = [Position.latitude, Position.longitude, Position.altitude];
lla0 = lla(1,:);

xyzNED = lla2ned(lla,lla0,'flat');