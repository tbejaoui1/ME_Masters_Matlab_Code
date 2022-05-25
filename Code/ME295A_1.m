%ME295A, Phase 1 Project
%Student Name: Takoua Bejaoui
%ID#:007576192
%Due Date: 12/18/21

%Collected static noise of the sensors via Matlab Mobile App. 
%This information is useful to determine whether or not the sensors' noise
%have a gaussian distribution or not. 
%Also gives you an idea on the variances of the noise as well. 
%Used iPhone8.
%iPhone seated on table with screen facing the ceiling. 
%I previously lost the previous script where I did everything from the
%original raw data, so I'm just using the transformed xcel file in this
%reuploaded version. 
%Saved each sensor's data in seperate .mat files before loaded into cloud.

%% Position
%load static position sensor log data from cloud
ld_pn = load('sensorlog_pos.mat'); 

%Position is the header name that the GPS outputs, while lat/long/alt are
%its features in their respective order.
Pos_nos= [ld_pn.Position.latitude, ld_pn.Position.longitude, ld_pn.Position.altitude, ld_pn.Position.hacc];

%Writing said table/matrix to an excel file
writematrix(Pos_nos, 'Pos_static_noise_ii.xlsx');
%Pos_nos = readtable('Pos_static_noise.xlsx');

t_sp = [1:length(Pos_nos)];

figure(1);
plot(t_sp,Pos_nos(:,1),t_sp, Pos_nos(:,2), t_sp, Pos_nos(:,3), t_sp, Pos_nos(:,4))
title('Static Noise of GPS')
xlabel('time [s]')
ylabel('decimal degrees (dd), meters (m)')
legend('Latitude noise','Longitude noise','altitude noise', 'Horizontal accuracy noise')

figure(2);
histogram(Pos_nos(:,1),'FaceColor','r')
title('Static Noise of GPS')
xlabel('Latitude values (dd)')
ylabel('frequency')

figure(3);
histogram(Pos_nos(:,2),'FaceColor','g')
title('Static Noise of GPS')
xlabel('Longitude values (dd)')
ylabel('frequency')

figure(4);
histogram(Pos_nos(:,3),'FaceColor','b')
title('Static Noise of GPS')
xlabel('Altitude values (m)')
ylabel('frequency')

figure(5);
histogram(Pos_nos(:,4),'FaceColor','m')
title('Static Noise of GPS')
xlabel('Horizontal Accuracy values (m)')
ylabel('frequency')
%% Accelerometer
%load static accelerometer data
%load static accelerometer sensor log data from cloud
ld_An = load('sensorlog_accel.mat'); 

%Accleration is the header name that the accelerometer sensor outputs, while X/Y/Z are
%its features in their respective order.
Accel_nos= [ld_An.Acceleration.X, ld_An.Acceleration.Y, ld_An.Acceleration.Z];

%Writing said table/matrix to an excel file
writematrix(Accel_nos, 'Accel_static_noise_ii.xlsx');


t_sa = [1:length(Accel_nos)];

figure(6);
plot(t_sa,Accel_nos(:,1),t_sa, Accel_nos(:,2), t_sa, Accel_nos(:,3))
title('Static Noise of Accelerometer')
xlabel('time [s]')
ylabel('m/s^2')
legend('X noise','Y noise','Z noise')

figure(7);
histogram(Accel_nos(:,1),'FaceColor','r')
title('Static Noise of Accelerometer')
xlabel('Accel X values (m/s^2)')
ylabel('frequency')

figure(8);
histogram(Accel_nos(:,2),'FaceColor','g')
title('Static Noise of Accelerometer')
xlabel('Accel Y values (m/s^2)')
ylabel('frequency')

figure(9);
histogram(Accel_nos(:,3),'FaceColor','b')
title('Static Noise of Accelerometer')
xlabel('Accel Z values (m^s^2)')
ylabel('frequency')
%% Magnetometer
%load static magnetometer data

%load static magnetometer sensor log data from cloud
ld_Mn = load('sensorlog_mag.mat'); 

%MagneticField is the header name that the magnetometer sensor outputs, while X/Y/Z are
%its features in their respective order.
Mag_nos= [ld_Mn.MagneticField.X, ld_Mn.MagneticField.Y, ld_Mn.MagneticField.Z];

%Writing said table/matrix to an excel file
writematrix(Mag_nos, 'Mag_static_noise_ii.xlsx');


t_sm = [1:length(Mag_nos)];

figure(10);
plot(t_sm,Mag_nos(:,1),t_sm, Mag_nos(:,2), t_sm, Mag_nos(:,3))
title('Static Noise of Magnetometer')
xlabel('time [s]')
ylabel('uT')
legend('X noise','Y noise','Z noise')

figure(11);
histogram(Mag_nos(:,1),'FaceColor','b')
title('Static Noise of Magnetometer')
xlabel('Magnetometer X values (uT)')
ylabel('frequency')

figure(12);
histogram(Mag_nos(:,2),'FaceColor','r')
title('Static Noise of Magnetometer')
xlabel('Magnetometer Y values (uT)')
ylabel('frequency')

figure(13);
histogram(Mag_nos(:,3),'FaceColor','g')
title('Static Noise of Magnetometer')
xlabel('Magnetometer Z values (uT)')
ylabel('frequency')

%% Gyro 
%load static gyro data

%load static gyroscope sensor log data from cloud
ld_Gn = load('sensorlog_angvel.mat'); 

%Angular velocity is the header name that the gyro sensor outputs, while X/Y/Z are
%its features in their respective order.
Gyro_nos= [ld_Gn.AngularVelocity.X, ld_Gn.AngularVelocity.Y, ld_Gn.AngularVelocity.Z];

%Writing said table/matrix to an excel file
writematrix(Gyro_nos, 'Gyro_static_noise_ii.xlsx');


t_sg = [1:length(Gyro_nos)];

figure(14);
plot(t_sg,Gyro_nos(:,1),t_sg, Gyro_nos(:,2), t_sg, Gyro_nos(:,3))
title('Static Noise of Gyro')
xlabel('time [s]')
ylabel('rad/s')
legend('X noise','Y noise','Z noise')

figure(15);
histogram(Gyro_nos(:,1),'FaceColor','b')
title('Static Noise of Gyro')
xlabel('Gyro X values (rad/s)')
ylabel('frequency')

figure(16);
histogram(Gyro_nos(:,1),'FaceColor','r')
title('Static Noise of Gyro')
xlabel('Gyro Y values (rad/s)')
ylabel('frequency')

figure(17);
histogram(Gyro_nos(:,1),'FaceColor','g')
title('Static Noise of Gyro')
xlabel('Gyro Z values (rad/s)')
ylabel('frequency')


%% Orientation (mix of Gyro + Accelerometer)
%load static orientation data

%load static gyroscope sensor log data from cloud
ld_On = load('sensorlog_orient.mat'); 

%Orientation is the header name that the Orientation sensor outputs, while X/Y/Z are
%its features in their respective order.
%X = Azimuth
%Y = Pitch
%Z = Roll
Orient_nos= [ld_On.Orientation.X, ld_On.Orientation.Y, ld_On.Orientation.Z];

%Writing said table/matrix to an excel file
writematrix(Orient_nos, 'Orient_static_noise_ii.xlsx');


t_so = [1:length(Orient_nos)];

figure(18);
plot(t_so,Orient_nos(:,1),t_so, Orient_nos(:,2), t_so, Orient_nos(:,3))
title('Static Noise of Orientation')
xlabel('time [s]')
ylabel('degrees')
legend('Azimuth noise','Pitch noise','Roll noise')

figure(19);
histogram(Orient_nos(:,1),'FaceColor','r')
title('Static Noise of Orientation Sensor')
xlabel('Azimuth values (deg)')
ylabel('frequency')

figure(20);
histogram(Orient_nos(:,2),'FaceColor','g')
title('Static Noise of Orientation Sensor')
xlabel('Pitch values (deg)')
ylabel('frequency')

figure(21);
histogram(Orient_nos(:,3),'FaceColor','b')
title('Static Noise of Orientation Sensor')
xlabel('Roll values (deg)')
ylabel('frequency')
