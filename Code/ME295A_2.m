%ME295A Collecting Sensor Data (GPS + Accel + Gyro(Orient + Angul.Vel) + Mag + GPS)
%Name: Takoua Bejaoui
%Due Date: 11/26/21

%Problem3: Plotting all the critical parameters while walking around 1 block.
%phone is held horizontally (screen face up). 



%% Collecting All Sensor Data Gyro + Accel + Mag + GPS on Smart Phone In Motion
clear  all;  close  all;  clc 

%walking and collecting all 
%variables of interest are lat(x), long(y), azimuth (deg),  angular
%velocity about azimuth axis, accel(x), accel(y)

%Loads the timetable file below from the app's sensor logs folder
load('sensorlog_everythang.mat'); 

%converting timetable to a table

T_Position1= timetable2table(Position);

T_Orientation1= timetable2table(Orientation);

T_AngularVelocity1= timetable2table(AngularVelocity);

T_Acceleration1= timetable2table(Acceleration);

T_Magnetometer1 = timetable2table(MagneticField);

%however using datenum and then diff, gets the time difference btwn
%timestamps
Np = datenum(T_Position1.Timestamp);  %creates a serial # to each timestamp
dt_secp = diff(Np)*24*60*60    ;          %finds the diff and converts each diff in sec


%however using datenum and then diff, gets the time difference btwn
%timestamps
No = datenum(T_Orientation1.Timestamp);  %creates a serial # to each timestamp
dt_seco = diff(No)*24*60*60    ;          %finds the diff and converts each diff in sec


%however using datenum and then diff, gets the time difference btwn
%timestamps
Nav = datenum(T_AngularVelocity1.Timestamp);  %creates a serial # to each timestamp
dt_secav = diff(Nav)*24*60*60    ;          %finds the diff and converts each diff in sec

%however using datenum and then diff, gets the time difference btwn
%timestamps
Nac = datenum(T_Acceleration1.Timestamp);  %creates a serial # to each timestamp
dt_secac = diff(Nac)*24*60*60    ;          %finds the diff and converts each diff in sec

%however using datenum and then diff, gets the time difference btwn
%timestamps
Nmag = datenum(T_Magnetometer1.Timestamp);  %creates a serial # to each timestamp
dt_secmag = diff(Nmag)*24*60*60    ;          %finds the diff and converts each diff in sec

%need to add a time vector:
%however there are diff data points for each feature, thus must create a
%different length time vector for each critical parameter within their respective
%automated time tables.
t_pos1 = [0:1:663]';
t_orient1 = [0:0.1:659.0]';
t_angvel1 = [0:0.1:667.2]';
t_accel1 = [0:0.1:667.2]';
t_mag1 = [0:0.1:658.9]';


%add the time vector to each table
T_Position1.t_pos1 = t_pos1;

T_Orientation1.t_orient1 = t_orient1;

T_AngularVelocity1.t_angvel1 = t_angvel1;

T_Acceleration1.t_accel1 = t_accel1;

T_Magnetometer1.t_mag1 = t_mag1;


%Writing said table to an excel file
writetable(T_Position1, 'Position_mobile.xlsx');
writetable(T_Orientation1, 'Orientation_mobile.xlsx');
writetable(T_AngularVelocity1, 'AngularVelocity_mobile.xlsx');
writetable(T_Acceleration1, 'Acceleration_mobile.xlsx');
writetable(T_Magnetometer1, 'Magnetometer_mobile.xlsx');



%Reads and grabs the selected data points wrt cell location in xcel
%However, won't grab the 1st column if it's in a timestamp format
Pos_Mobile= xlsread('Position_mobile.xlsx','Sheet1','A1:H665');
Orient_Mobile= xlsread('Orientation_mobile.xlsx','Sheet1','A1:E6592');
AngVel_Mobile= xlsread('AngularVelocity_mobile.xlsx','Sheet1','A1:E6674');
Accel_Mobile= xlsread('Acceleration_mobile.xlsx','Sheet1','A1:E6674');
Mag_Mobile = xlsread('Magnetometer_mobile.xlsx', 'Sheet1', 'A1:E6591');

lat_pos= Pos_Mobile(:,1);        %  latitude
long_pos= Pos_Mobile(:,2);     %longitude
alti_pos= Pos_Mobile(:,3);  % altitude (height)
speed = Pos_Mobile(:,4);    % speed (m/s)
course = Pos_Mobile(:,5);   %course (deg)
hac_pos = Pos_Mobile(:,6);  % horizontal accuracy
t_pos = Pos_Mobile(:,7);     %time in [s]

x_orient=Orient_Mobile(:,1);        % psi angle (azimuth)
y_orient=Orient_Mobile(:,2);     %theta angle (pitch)
z_orient=Orient_Mobile(:,3);  % phi angle (roll)
t_orient= Orient_Mobile(:,4);     %time in [s]

%the speed at which you yaw, pitch or roll ur phone
pitch_angvel= AngVel_Mobile(:,1);        % angvel about x axis (of phone)cw=+
roll_angvel= AngVel_Mobile(:,2);     % angvel about y axis (of phone)
yaw_angvel= AngVel_Mobile(:,3);   %angvel about z axis (of phone) where cc = +
t_angvel= AngVel_Mobile(:,4);     %time in [s]

x_accel=Accel_Mobile(:,1);        %  accel along x axis
y_accel=Accel_Mobile(:,2);     %accel along y axis
z_accel=Accel_Mobile(:,3);  %accel along z axis
t_accel= Accel_Mobile(:,4);     %time in [s]

Mag_X = Mag_Mobile(:,1);     %magnetic field, uT, intensity along East/West
%when closer to the south pole, this value becomes more or less negative
%but when closer to the north pole, this value becomes more positive
%negative = magnetic west, positive = magnetic east

Mag_Y = Mag_Mobile(:,2);     %magnetic field, uT, ambient strength, up or down wrt Earth. 
%Closer to the North pole brings out a large negative Y field, while a
%large positive Y magnetic value means it's closer to the South pole.
%the y-axis points towards the south pole but away from the north pole

%but values are opposite for magnetic north....north pole is true north.
%so positive is actually pointing towards magnetic north and negative when
%pointing towards magnetic south.

%positive = magnetic north, negative = magnetic south

Mag_Z = Mag_Mobile(:,3);     %magnetic field, uT, Up or Down
%closer to the south pole, this value becomes less negative
%closer to the north pole, this value becomes more negative
%negative = pointing down
%however, the way the phone is tilted changes which axis tells where
%magnetic north is at....for ex, tilting the phone vertically and then to
%the side, makes the x-axis determine magnetic north's intensity instead of the 
%y -axis as it was originally. 

%Regardless, it is the overall magnetic intensity of the x-y-z values given
%by the magnetometer that will allow us to calculate declination (angle between
% magnetic north and true or geographic north) and inclination angles 
%(angle between magnetic intensity and magnetic north).


%t_every = t_pos;   %since this is the max amount of time steps in comparison



% figure()   %all vectors need the same length
% plot(t_every,lat_pos,t_every,long_pos,t_every,yaw_angvel, t_every, x_orient, t_every, x_accel, t_every, y_accel, t_every, z_accel)
% title('Time vs. Critical Parameters of Phone')
% xlabel('Time [s]')
% ylabel('Critical Parameters [(lat)deg, (long)deg, (angvel_azi) rad/s, (orient_azi) deg, (accel_x) m^2/s, (accel_y) m^2/s, (accel_z) m^2/s]')
% legend('latitude (deg)','longitude (deg)','angular velocity azimuth (rad/s)', 'orient_azi (deg)', 'accel_x (m^2/s)', 'accel_y (m^2/s)', 'accel_y (m^2/s)');

%Only selected critical parameters are focused on
figure(1)
plot(t_pos,lat_pos,t_pos,long_pos, t_pos, speed) 
title('Time vs. Mobile Position Parameters of Phone (GPS)')
xlabel('Time [s]')
ylabel('Critical Parameters [(lat)deg, (long)deg, (velocity)m/s]')
legend('latitude (deg)','longitude (deg)', 'speed (m/s)');

figure(2)
plot(t_angvel,yaw_angvel) 
title('Time vs. Angular Velocity of Phone')
xlabel('Time [s]')
ylabel('Critical Parameters [(angvel_azi) rad/s]')
legend('angular velocity azimuth (rad/s)');

figure(3)
plot(t_orient, x_orient) 
title('Time vs. Orientation Parameter of Phone')
xlabel('Time [s]')
ylabel('Critical Parameters [(orient-azimuth) deg]')
legend('orient - azimuth (deg)');

figure(4)
plot(t_accel, x_accel, t_accel, y_accel, t_accel, z_accel)
title('Time vs. Acceleration Parameters of Phone')
xlabel('Time [s]')
ylabel('Critical Parameters [(accel-x) m/s^2, (accel-y) m/s^2, (accel-z) m/s^2]')
legend('accel-x (m/s^2)', 'accel-y (m/s^2)', 'accel-z (m/s^2)');
