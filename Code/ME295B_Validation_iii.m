%ME295B, Phase 2 Project
%Student Name: Takoua Bejaoui
%ID#:007576192
%Due Date: 05/10/22
clear all; close all;
%This script is Approach 2 for the discontinuous kalman observer model 
% for ME295B. And it is served to validate other runs besides the main one
% for the project.
%Where measured output is GPS in either the X and Y directions and the
%control inputs are acceleration in both the x and y directions. 
%It's advisable to run section by section.
%This is the master script for Approach 2, where absolute acceleration is
% captured via rotation matrix. Thus, gravitational component is derived
% from the rotational matrix at each data step. This is an attempt to 
% seperate true/linear acceleration in absolute reference frame. This is 
% to reduce the error offset of the kalman/observer model's estimated, 
% xhat, output.

%% Converting GPS's lla (long, lat and altitude) to NED/ENU

%Phone's GPS outputs latitude,longitude, altitude and speed.(llav) 
%where latitude = Y, longitude = X, altitude = a (or 'h' in literature), and speed = v

%There are two methods to accessing the raw .mat data; either converting
%the given timetable (default) to a table that can be accessed as an excel
%sheet; or access the given timetable as is. 
%If wanting to run code using the table, then I suggest using the script
% titled: ME295A_2 setup and using 'readtable' instead of 'load'. 
%The ME295A_2 contains the tables and can use the 'dot functions' to
%access the different classes/features/columns/headers within that table.

%But if running the given timetable of raw .mat data, then continue with
%the following setup. 


%loading the timetables for each type of parameter collected via the Matlab
%mobile app. The sensors will be listed as follows...
ld_m = load('accelwigglenewnorth1.mat'); 

%We only care for latitude (phi), longitude(lambda) and altitude(h)
lla = [ld_m.Position.latitude, ld_m.Position.longitude, ld_m.Position.altitude];
lla0 = lla(1,:);  %starting point in position mobile data - collected in file: ME295A_2
                   %also called P_ref (reference position) in literature
%used Matlab online sjsu (student) version to access necessary toolboxes:
%Navigation toolbox
%Sensor Fusion and Tracking toolbox
%UAV toolbox
%Have the option of using either ENU or NED
%xyzNED = lla2ned(lla,lla0,'flat');

%Will use ENU for phone
xyzENU = lla2enu(lla,lla0,'flat');

figure(1);
grid on
grid minor
%longitude along X and latitude along Y.
plot(lla(:,2), lla(:,1));

figure(2);
grid on
grid minor
%plots xyz data wrt given lat and long in ENU frame (x-East, y-North,
%z-Down)
plot(xyzENU(:,1), xyzENU(:,2))

%% Getting Yaw Readings thru Accelerometer/Magnetometer via 'ecompass'

%magnetometer readings are strength values pertaining to how close or how
%far one is from either the magnetic North or South poles. 
%For instance, if the reference frame is ENU for the phone, then a higher
%positive Y magnetic value (tesla units, uT) indicates that you're in the
%N-E or N-W quadrants; while a negative Y magnetic value indicates that you're in the
% S-E and S-W quadrants. If X is at 0 uT, then it is pointing exactly East or West and Y
%should be pointing towards magnetic north or magnetic South, respectively.
%A Y = 0uT is also pointing at East or West. 
%if phone is tilted so that X axis is pointing towards the North, it reads
%a high positive value. 
%If phone is tilted so that X axis is pointing towards the South, it reads
%a high negative value. 
%It seems like X and Y take turns in determining where magnetic North and South
%poles are wrt orientation of phone. (basically where the X and Y axes of
%the body frame/phone is pointing at.)
% Z, in the magenetometer is always pointing down (not south) normal to the
% x-y plane. And is usually a negative value, since the reference frame is
% ENU for the phone. And its average value is usually around the value of
% the latitude that ur on .

%For the accelerometer data, due to sensor noise and the way the phone is
%horizonatally unstably held throughout the data collection process, the z
%axis values fluctuate about the true 9.8m/s^2 value
%Y value is the acceleration of motion of the device (basically myself as
%I'm holding the phone horizontally with the Y-axis pointing forward.)

%loading the magnetometer and accelerometer data.
%proper accelerometer is accelerometer data

%Listing the data based on values and headers
Allmagreads = [ld_m.MagneticField.X, ld_m.MagneticField.Y, ld_m.MagneticField.Z];
AllAccelreads = [ld_m.Acceleration.X, ld_m.Acceleration.Y, ld_m.Acceleration.Z];

Mag_x = Allmagreads(:,1);
Mag_y = Allmagreads(:,2);
Mag_z = Allmagreads(:,3);

Accel_x = AllAccelreads(:,1);
Accel_y = AllAccelreads(:,2);
Accel_z = AllAccelreads(:,3);


%In order to use the ecompass library, the time frames, btwn
% the accelerometer and magnetometer need to be equal. Thus, interpolation
%is used. 

%from the initial time to its end wrt the accelerometer data (since
%it has more data points than the magnetometer).
%Looking at the timestamps of the raw data helped determine the time 
%interval. 
%Using just one of the variables is fine, since Y, Z and X are all the same
%size
t_Accel = [0:0.1:(length(Accel_x)/10)-1]';

% x2-x1/n-1, where linspace(x1, x2, n) generates n points that are evenly
% spaced.
t_Magn = linspace(0,t_Accel(end), length(Mag_x))';


%Linear Interpolation of the Accelerometer and Magnetometer signals
Magnetometer_inter_x = interp1(t_Magn, Mag_x,t_Accel);
Magnetometer_inter_y = interp1(t_Magn, Mag_y,t_Accel);
Magnetometer_inter_z = interp1(t_Magn, Mag_z,t_Accel);

%Have to do similar to the accelerometer itself...
t_accel_raw = linspace(0,t_Accel(end), length(Accel_x))';
Accelerometer_inter_x = interp1(t_accel_raw, Accel_x,t_Accel);
Accelerometer_inter_y = interp1(t_accel_raw, Accel_y,t_Accel);
Accelerometer_inter_z = interp1(t_accel_raw, Accel_z,t_Accel);

Allmagreads_inter = [Magnetometer_inter_x, Magnetometer_inter_y, Magnetometer_inter_z];
AllAccelreads_inter = [Accelerometer_inter_x,Accelerometer_inter_y,Accelerometer_inter_z];

magneticFieldStrength = Allmagreads_inter(:,:);  %nth reading from magnetometer
Acceleration = AllAccelreads_inter(:,:);         %nth reading from the accelerometer

%then using the ecompass formula:
%which gives the rotation information, either wrt to a quaternion or a
%rotation matrix, from the inertial local ENU/NED (in this case, ENU) to
%the child frame (sensor/body frame) 
% This formula will give you declination if the device is pointed towards
% true north. This is b/c the magnetometer only senses magnetic north and
% south poles. Thus the difference between magnetic north (MN) and true
% north (TN) will be the declination angle. 
%But pointing the device anywhere else will just give me its position in
%relation to magnetic poles. 
%Its output would be a rotation vector, wrt each sample reading, either in
%quaternions or via a rotation matrix. 
orientation = ecompass(Acceleration,magneticFieldStrength,'quaternion','ReferenceFrame','ENU');
%orientation1 = ecompass(Acce6290,magneticFieldStrength6290,'rotmat','ReferenceFrame','ENU')

%% Yielding Yaw angle wrt Magnetic North via Eulerd Class

%convert quaternion orientation/rotation vector to Euler rotation angles:
%ZXY is selected b/c according to my sensor's reference frame:
%rotation about Z will give yaw angle, rotation about X will give me pitch
%and rotation about Y will give me roll
%the order of this is important b/c matlab's orientation, as well as most
%outputs pertaining to the fused orientation will yield: [yaw, pitch, roll]
eulerAnglesDegrees = eulerd(orientation,'ZXY','frame');
%therefore yaw = eulerAnglesDegree(1) 
%yaw = eulerAnglesDegrees(1)   %The first output is yaw
%roll = eulerAnglesDegrees(3)  %The last output is a roll

%Another way of getting yaw, pitch and roll from quaternion:
%But only takes in the coefficients of the real and imaginary components
%so, for ex: orientation = [0.1370   -0.8494   -0.5067    0.0544]
%orientation2 = [0.10225,-0.9890,0.0643,-0.0856] (some point:6290.)
%[yaww, pitch, roll] = quat2angle(coeff_quat, "ZXY")
%yawd = rad2deg(yaww)  %it's the same as the values for eulerAnglesDegrees

%Assuming YAW = eulerAnglesDegrees(1)
%At datapoint 1632 --> 171.6 deg
%At datapoint 2715 --> -40.5 deg
%At datapoint 4399 --> 8.85 deg
%At datapoint 5633 --> 49.37 deg

%angles without the rotational fixes as is depicted in eulerAnglesDeg.

%This method treats N = 0deg, E = 90 deg, S=+/-180 deg, W = -90deg


%% Plotting the Yaw angles of 'ecompass' EulerAngleDeg variable wrt True North
%To convert yaw wrt magnetic north to true north, must consider the
%declination angle. 
%13 deg 3 min is the declination wrt where the data was geographically
% collected. It is subtracted due to being + East.
%Thus, 3 min/60 = 05deg.
%Total declination in deg: 13.05 deg
psi_d = eulerAnglesDegrees(:,1) - 13.05;

figure(3);
plot(psi_d)
xlabel('Time Steps (at every 0.1 sec)')
ylabel('Psi angle wrt True North [deg]')


%% Processing Accelerometer, Magnetometer and GPS data before Kalman Filter input

%Accelerometer
a_x= Accel_x;
a_y = Accel_y;

t_a = [0:0.1:(length(a_x)/10)-1]';

%GPS
X_GPS_raw = xyzENU(:,1);
Y_GPS_raw = xyzENU(:,2);

t_GPS = linspace(0,t_a(end), length(X_GPS_raw))';

%Yaw angle wrt True North
psi_r_raw = psi_d*(pi/180);

t_psi = linspace(0,t_a(end), length(psi_r_raw))';

%Linear Interpolation of the GPS and Magnetometer signals

t= t_a;

X_GPS = interp1(t_GPS, X_GPS_raw,t);
Y_GPS = interp1(t_GPS, Y_GPS_raw,t);
psi_r = interp1(t_psi, psi_r_raw,t);

figure(4)
plot(X_GPS, Y_GPS)

%check if it matches with the xyz_ENU_x and xyz_ENU_y plotted previously.

%% Setting up the State Matrices for dqle (LQR estimated descritization) Method

%[L, P, Z, E] = dqle(A,G,C,Q,R);

%discrete (discontinous) linear quadratic estimater for Kalman Filter
% requires the following matrices: (A, G, C, Q, R).
%Please check report for mathematical derivations wrt FBD...
%Where A = state transition matrix of the coefficient without the inputs
%G = coefficient matrix for the input (accelerometer) noise
%C = measured/observable selected state's(s') coefficient matrix
%Q = process noise variance matrix wrt states
%R = output measured noise variance matrix (GPS)

%With respect to the report:
%With respect to X_GPS as the measured output and X_hat_k, X_hat_k+1 as
% the estimated outputs, which are determined by raw accelerometer data:
dt = 0.1; %thnx to the interpolating method, both output (GPS) and
%input (accelerometer) have the same time intervals. 
A = [0,1;-1,2];
%overall control coefficent matrix along the x axis:
Bx = [0,0;(dt^2).*cos(psi_r), (dt^2).*sin(psi_r)];
%overall control coefficient matrix along the y axis:
By = [0,0;-(dt^2).*sin(psi_r), (dt^2).*cos(psi_r)];

C = [0,1];
Q = 0.15; %(variance of the process noise wrt 
% accelerometer's static data(when in motion), where accelerometers -x and -y have been 
% equated to have around equal variances.) Thus (var(a_x ==var(a_y)).

R = 4; %standard deviation/accuracy is +/- 4.9m. Observable/measured 
%state's noise's variance (GPS in either the X or Y axis). R can be in a
%range from 1 : 25. 1 and 4 are most accurate, b/c the resolution of the
%GPS sensor used was at every 1m or 2 m. Witnessed this in the model output
%of xhat. 
G = [0;dt^2]; %noise coefficient matrix wrt accelerometer's process
%noise. The addition or subtraction of cos and sin terms yields a periodic
%curve  - in this case a sine wave with a slight bias. Due to the methods
%used where G is considered to be a 'static' matrix, this noise coefficient
%matrix has been simplified to this form. Thus removing the non linear sin
%and cos terms. The static noise histogram plots captured in ME295A_1
%script indicates that there is indeed a bias in the accelerometer
%static/stationary data (where mean, for both x and y, are not equal to
% zero). Thus, acceleration noise is just grouped as one entity of noise,
% rather than being distinguished as n_x and n_y along the accelerometer's
% x and y axis, respectively. 

%Remember, Q and R are variance matrices that measure 'confidence'
% not 'accuracy' of the sensor input and desired measured output,
% respectively. 
%To get the Q matrix, need to grab the variance of the static noise data
%from the phone's accelerometer. Thus, as the phone is statically placed
%horizontally on the table, log stream the data for a few minutes, and then
%calculate the variance of that. This was done in the ME_295A_1 script.

%Calculated the R values in ME_295A_1
%Didn't have to b/c it would lead to an overconfident matrix due
%to the GPS's ability 
%R_x = [3.136363636126141e-11];
%R_y = [3.379113636223254e-09];

%We'll compare the R approaches....
%The R_x and R_y yield Eigenvalues closer to zero, but still unstable.
%Highter R values leads to values closer to + 1 but still have oscillatory
%behavior due to the imaginary factor. 
%It has been advised not to use a computed GPS's noise variance since it
%will yield and over confident matrix. Thus, R is based on the overall
%accuracy as stated per manufacturer's specs. 

%Anyways, using the dqle (it's been changed to dlqe)

[L,P,Z,E]=dlqe(A,G,C,Q,R);
% Where L = is the observer gain (like the 'K' gain in feedback loop)
%P = variance values of state (Ricatti soln)
%Z = variance values of measured (error covariance)
%E = eigenvalues (poles), which gives an idea of the stability of the
%system.

%E are two unstable values with some oscillations. positive real +
%imaginary parts that are close to zero.

% P is a 2x2 matrix, with positive values close to zero. Thus, P is stable.

%Z is a 2x2 matrix with positive values close to zero, as well. 

%Testing the eigenvalues of the A matrix
eig_A = eig(A); %two positive poles, unstable open system; however, this is
%not a continous system. A discontinuous system behaves according to the
%unit circle domain. Thus as long as the values are within the unit circle,
% -1:1, then system is stable.

%Since the Kalman filter is not a controller but an observer does it still
%make sense to test for controllability and observability of the system?

%Cont_x = ctrb(A,Bx); couldn't do it due to the size of the B matrix
%The B matrix isn't constant.

%Let's test for observability
Ob = obsv(A,C);
rank_Ob = rank(Ob); %Must equal number of rows to be full rank, 2
%system is fully observable.

%Since the eigenvalues, E, are with respect to a discontinous kalman
%filter, this system is stable. All eigenvalues must be between values of
%1, -1. Due to the unit circle domain. 

 %Modification of Q and R matrices really effects the final L value.
 %Ideally want a different L (observer gain) value for x and y axes. 

%% Running the Discontinuous Observer Model/Kalman Filter
%With respect to the X parameter of the model.(xhats,X_GPS).

%Using the interpolated accelerometer values as inputs for filter
a_xx = Accelerometer_inter_x;
a_yy = Accelerometer_inter_y;
a_zz = Accelerometer_inter_z;

%stating initial conditions for all states and estimated states
x_hat = [0;0]; %initial conditions for all the estimated output states (wrt X).
           %wrt X_hat(k-1) and X_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);
for k = 1:k_end-1
    %solving for Bx in the for loop instead
    B_x = [0,0;(dt^2).*cos(psi_r(k)), (dt^2).*sin(psi_r(k))];
    u = [a_xx(k,1); a_yy(k,1)];

    %Measured X_GPS raw
    z = X_GPS(k,1);

    %one of the output states wrt the previous one
    %the generic full kalman filter general equation as stated in 
    %the report
    %Note: multiplying L by 0 gives an estimated output, xhat, that
    %experiences a large offshoot. This is to be expected since the
    %accelerometer's double integration, which yeilds position, is noisy.
    x_hat(:,k+1) = A*x_hat(:,k) + B_x*u+L*(z-C*x_hat(:,k));
    
end

%estimated outputs plotted, x(k-1) and x(k)
figure(5); plot([1:k_end],x_hat)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('X')
legend('x(k-1)', 'x(k)')
%the earlier curve is x(k-1) and the latter is x(k)

%plotted raw GPS
figure(6); plot([1:k_end],X_GPS)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('X')
legend('raw X-GPS')


%plotted both curves to compare
figure(7); plot([1:k_end], x_hat, [1:k_end], X_GPS)
xlabel('Time Steps (at every 0.1 sec)'); 
ylabel('Comparison btwn Estim. and Meas. X')
legend('x(k-1)', 'x(k)','raw X-GPS')



%% Running the Discontinuous Observer Model/Kalman Filter wrt absolute accel only
%with respect to the X and Y parameters of the model....

%attempt to get absolute acceleration along the z axis along the ENU frame.
% The linear acceleration is captured then, by removing gravitational component.
%Latter is done via multiplying the transpose of the rotation matrix to the raw
%interpolated accelerometer data. 
%tilt angles occur due to a sudden roll (x-z plane) or pitch (y-z plane)
%of the device, in this case phone. Since phone was only laid flat
%horizontally on palm, unexpected tilts occurred. 

% [Accel_E, Accel_N, Accel_U] = [R.T][acc_x, acc_y, acc_z]

orientation1 = ecompass(Acceleration,magneticFieldStrength,'rotmat','ReferenceFrame','ENU');

%Test with one, 3x3, rotation matrix for the first acceleration data sample
%absol_accel_ENU1 = orientation1(:,:,1)'*[a_xx(1,1);a_yy(1,1); a_zz(1,1)];

%Now doing for all the values
i_end = length(t_a);
for i = 1: i_end
    absol_accel_ENU(:,1,i) = orientation1(:,:,i)'*[a_xx(i,1);a_yy(i,1); a_zz(i,1)];
end

%The for loop above showed the gravitational component along all axes.
%Thus, only the z-axis had a non zero value. 
%This value needs to be subtracted from the tilt angle formulations derived
%in the ME295B_F_A1, instead of just using 9.81 as the g component.

%done with respect to y axis
%pitch tilt angle
theta_tilt_y2 = atan2(a_zz,a_yy); 
% figure; plot(theta_tilt*180/pi)
z1_end = length(t_a);

for z1 = 1:z1_end
    z_g = absol_accel_ENU(3,:,:); % 1x1x6664
    a_yy2(z1) = a_yy(z1) + z_g(1,1,z1).*cos(theta_tilt_y2(z1,1)); 
end
%figure(8); plot(a_yy2)
%the plot shows acceleration from zero axis and below.


%done with respect to acceleration in the y direction
%roll tilt angle
theta_tilt_x2 = atan2(a_zz,a_xx); 
for z1 = 1:z1_end
    z_g = absol_accel_ENU(3,:,:); % 1x1x6664
    a_xx2(z1) = a_xx(z1) + z_g(1,1,z1).*cos(theta_tilt_x2(z1,1)); 
end
%figure(9); plot(a_xx2)
%the plot indicates that acceleration is about the zero axis


%% Running the Discontinuous Observer Model/Kalman Filter wrt linear accel only - Continue

%stating initial conditions for all states and estimated states
x_hat = [0;0]; %initial conditions for all the estimated output states (wrt X).
           %wrt X_hat(k-1) and X_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);
for k = 1:k_end-1
    %solving for Bx in the for loop instead
    B_x = [0,0;(dt^2).*cos(psi_r(k)), (dt^2).*sin(psi_r(k))];
    a_xx3 = a_xx2';
    a_yy3 = a_yy2';
    u = [a_xx3(k,1); a_yy3(k,1)];

    %Measured X_GPS raw
    z = X_GPS(k,1);

    %one of the output states wrt the previous one
    %the generic full kalman filter general equation as stated in 
    %the report
    %Note: multiplying L by 0 gives an estimated output, xhat, that
    %experiences a large offshoot. This is to be expected since the
    %accelerometer's double integration, which yeilds position, is noisy.
    x_hat(:,k+1) = A*x_hat(:,k) + B_x*u+L*(z-C*x_hat(:,k));
    
end

%estimated outputs plotted, x(k-1) and x(k)
figure(10); plot([1:k_end],x_hat)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('X')
legend('x(k-1)', 'x(k)')

%plotted raw GPS
figure(11); plot([1:k_end],X_GPS)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('X')
legend('raw X-GPS')

%plotted both curves to compare
figure(12); plot([1:k_end], x_hat, [1:k_end], X_GPS)
xlabel('Time Steps (at every 0.1 sec)'); 
ylabel('Comparison btwn Estim. and Meas. X')
legend('x(k-1)', 'x(k)','raw X-GPS')

%% Running the Discontinuous Observer Model/Kalman Filter
%With respect to the Y parameters of the model.(yhats,Y_GPS).

%Using the interpolated accelerometer values as inputs for filter
a_xx = Accelerometer_inter_x;
a_yy = Accelerometer_inter_y;
a_zz = Accelerometer_inter_z;



%stating initial conditions for all states and estimated states
y_hat = [0;0]; %initial conditions for all the estimated output states (wrt Y).
           %wrt Y_hat(k-1) and Y_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);
for k = 1:k_end-1
    %solving for Bx in the for loop instead
    B_y = [0,0;-(dt^2).*sin(psi_r(k)), (dt^2).*cos(psi_r(k))];

    u = [a_xx(k,1); a_yy(k,1)];

    %Measured X_GPS raw
    z = Y_GPS(k,1);

    %one of the output states wrt the previous one
    %the generic full kalman filter general equation as stated in 
    %the report
    %Note: multiplying L by 0 gives an estimated output, xhat, that
    %experiences a large offshoot. This is to be expected since the
    %accelerometer's double integration, which yeilds position, is noisy.
    y_hat(:,k+1) = A*y_hat(:,k) + B_y*u+L*(z-C*y_hat(:,k));
    
end

%estimated outputs plotted, y(k-1) and y(k)
figure(13); plot([1:k_end],y_hat)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('Y')
legend('y(k-1)', 'y(k)','location', 'north')
%the earlier curve is y(k-1) and the latter is y(k)

%plotted raw GPS
figure(14); plot([1:k_end],Y_GPS)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('Y')
legend('raw Y-GPS','location', 'north')

%plotted both curves to compare
figure(15); plot([1:k_end], y_hat, [1:k_end], Y_GPS)
xlabel('Time Steps (at every 0.1 sec)'); 
ylabel('Comparison btwn Estim. and Meas. Y')
legend('y(k-1)', 'y(k)','raw Y-GPS','location', 'north')



%% Running the Discontinuous Observer Model/Kalman Filter wrt linear accel only - Continue
%with respect to the Y parameters of the model....

%stating initial conditions for all states and estimated states
y_hat = [0;0]; %initial conditions for all the estimated output states (wrt Y).
           %wrt Y_hat(k-1) and Y_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);

for k = 1:k_end-1
    %solving for Bx in the for loop instead
    B_y = [0,0;-(dt^2).*sin(psi_r(k)), (dt^2).*cos(psi_r(k))];
    a_xx3 = a_xx2';
    a_yy3 = a_yy2';
    u = [a_xx3(k,1); a_yy3(k,1)];

    %Measured X_GPS raw
    z = Y_GPS(k,1);

    %one of the output states wrt the previous one
    %the generic full kalman filter general equation as stated in 
    %the report
    %Note: multiplying L by 0 gives an estimated output, xhat, that
    %experiences a large offshoot. This is to be expected since the
    %accelerometer's double integration, which yeilds position, is noisy.
    y_hat(:,k+1) = A*y_hat(:,k) + B_y*u+L*(z-C*y_hat(:,k));
    
end

%estimated outputs plotted, y(k-1) and y(k)
figure(16); plot([1:k_end],y_hat)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('Y')
legend('y(k-1)', 'y(k)','location', 'north')

%plotted raw GPS
figure(17); plot([1:k_end],Y_GPS)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('Y')
legend('raw Y-GPS','location', 'north')

%plotted both curves to compare
figure(18); plot([1:k_end], x_hat, [1:k_end], X_GPS)
xlabel('Time Steps (at every 0.1 sec)');
ylabel('Comparison btwn Estim. and Meas. X')
legend('y(k-1)', 'y(k)','raw Y-GPS','location', 'north')

%% Final Plotting of X vs Y

figure(19);
plot(x_hat, y_hat, X_GPS, Y_GPS)
xlabel('X'); ylabel('Y')
legend('y-Kalman vs.x-Kalman','Y-GPS vs. X-GPS', 'location', 'southwest' )

%make sure the aspect ratio is true wrt data. 
figure(20);
plot(x_hat(2,:), y_hat(2,:),'.', X_GPS, Y_GPS,'linewidth',2)
xlabel('X'); ylabel('Y')
%pbaspect([1 1 1])
daspect([1 1 1])
xlim([-10,10])
legend('y-Kalman vs x-Kalman','Y-GPS vs X-GPS', 'location', 'southwest' )

%% Low Pass Filter Comparisons
%Along the X axis.
%alpha factor to be tuned.
%stating initial conditions for all states and estimated states
x_filt = 0; %initial conditions for all the estimated output states (wrt X).
           %wrt X_hat(k-1) and X_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);
alpha = 0.95;
for k = 1:k_end-1
    x_filt(k+1) = alpha*x_filt(k) + (1-alpha)*X_GPS(k);
    
end

%estimated outputs plotted,x(k) only
figure(21); plot([1:k_end],x_hat(2,:), [1:k_end],x_filt, [1:k_end],X_GPS )
xlabel('Time Steps (at every 0.1 sec)'); ylabel('X')
legend('x-Kalman','x-lpf', 'X-GPS', 'north' ) %lpf = low pass filter

%% Low Pass Filter Comparisons
%Along the Y axis.
%alpha factor to be tuned.
%stating initial conditions for all states and estimated states
y_filt = 0; %initial conditions for all the estimated output states (wrt X).
           %wrt Y_hat(k-1) and Y_hat(k)

%y0 = [0;0]; %initial conditions for all the states (wrt Y).
%yhat0= [0;0]; %initial cond, for state observer block diagram's integrator
k_end = length(t_a);
alpha = 0.95;
for k = 1:k_end-1
    y_filt(k+1) = alpha*y_filt(k) + (1-alpha)*Y_GPS(k);
    
end

%estimated outputs plotted,x(k) only
figure(22); plot([1:k_end],y_hat(2,:), [1:k_end],y_filt, [1:k_end],Y_GPS )
xlabel('Time Steps (at every 0.1 sec)'); ylabel('Y')
legend('y-Kalman','y-lpf', 'Y-GPS', 'location','southeast' ) %lpf = low pass filter

%% Low Pass Filter Comparisons
% x and y

figure(23); plot(x_hat(2,:),y_hat(2,:), x_filt,y_filt, X_GPS,Y_GPS )
xlabel('X'); ylabel('Y')
legend('y-Kalman vs.x-Kalman','x-lpf vs. y-lpf', 'X-GPS vs. Y-GPS', 'location','southwest' ) %lpf = low pass filter
daspect([1 1 1])
xlim([-10,10])
%% Error Curve Comparisons in either axes.

%Plotting error where error = reference - estimated feedback value
k_end = length(t_a);
for k = 1:k_end
    %
    e_x(k,1) = X_GPS(k,1) - x_hat(2,k);
    e_xf(k,1) = X_GPS(k,1) - x_filt(k);

    e_y(k,1) = Y_GPS(k,1) - y_hat(2,k);
    e_yf(k,1) = Y_GPS(k,1) - y_filt(k);
end

figure(24); plot([1:k_end],e_x,[1:k_end],e_xf)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('error along the x axis')
legend('error x-Kalman','error x-lpf') %lpf = low pass filter

figure(25); plot([1:k_end],e_y,[1:k_end], e_yf)
xlabel('Time Steps (at every 0.1 sec)'); ylabel('error along the y axis')
legend('error y-Kalman','error y-lpf') %lpf = low pass filter
