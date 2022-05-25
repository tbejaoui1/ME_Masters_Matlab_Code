%ME295B, Phase 2 Project
%Student Name: Takoua Bejaoui
%ID#:007576192
%Due Date: 05/10/22

%It's advisable to run section by section.
%This script is to run tests for capturing and validating the inputs for
%the control systems design. 
%Thus, sections that have 'Test' as part of their title weren't used in the
%final Master script.
%% State Space (SS) Representation of Sensor Fusion of GPS + Accelerometer 
%Sensor Fusion = Kalman Filter (is an observer not a controller)

%The data collected via Matlab Mobile App.
%Phone was positioned horizontally upon palm of hand. With Y-axis pointing
%forward, X-axis pointing East and Z-axis pointing upward.
%The lat and long of GPS can be converted from dd (lla format) to ECEF or
%to NED or ENU. (selected NED)
%All collected data wrt GPS and IMU can be found under MobileSensorData
%under file path: C:/Users/ilovehorses/Matlab Drive/MobileSensorData.
%Have two options to convert geodetic lla GPS raw data to local NED frame:
%either using Matlab's llatoNED function (need tools, but also can be done
%in Matlab Online with student account.
%Or:
%A series of hand calculations were done to convert geodetic lla to NED.
%The 2D representation has state variables : [x, y, psi], where x and y are
%long and lat respectively in m, and psi is the yaw angle of
%orientation in deg. 
%The local NED body frame converts the lla data of GPS to have its origina
%at the center of the robot,
%along with the changed local coordinate ref system wrt acceleration and
%orientation.
%Hand written work can be found in report as well as below.

%% Converting GPS's lla (long, lat and altitude) to NED/ENU
%Phone's GPS outputs latitude,longitude, altitude and speed.(llav) 
%where latitude = Y, longitude = X, altitude = a (or 'h' in literature), and speed = v

%run the script titled: ME295A_2 first before running this section. 
%b/c that script contains the tables and can use the 'dot functions' to
%access the different classes/features/columns/headers within that table.

%loading the table of the mobile position readings - collected from the previous script
%(ME295A_2)
Position = readtable('Position_mobile.xlsx');

%We only care for latitude (phi), longitude(lambda) and altitude(h)
lla = [Position.latitude, Position.longitude, Position.altitude];
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
Magnetometer = readtable('Magnetometer_mobile.xlsx');
Accelerometer = readtable('Acceleration_mobile.xlsx');

%Listing the data based on values and headers
Allmagreads = [Magnetometer.X, Magnetometer.Y, Magnetometer.Z];
AllAccelreads = [Accelerometer.X, Accelerometer.Y, Accelerometer.Z];

Mag_x = Allmagreads(:,1);
Mag_y = Allmagreads(:,2);
Mag_z = Allmagreads(:,3);

Accel_x = AllAccelreads(:,1);
Accel_y = AllAccelreads(:,2);
Accel_z = AllAccelreads(:,3);


%collecting the nth readings from each
%1632 data point (Y pointing S)
%2715 data point (Y pointing W)
%4399 data point (Y is pointing back towards N)
%5633 data point (Y is pointing E)

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
t_accel_raw = [Accelerometer.t_accel1];
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

%% ProcessPhoneData2 function - Test
% the purpose of this section is to test a function for the purposes
%of quickly interpolating sensor data and using it later with the 
%qfix*qtrue mathematics to get the corrected euler angles (esp. the Yaw)

%This way then it is possible to extract linear acceleration or to
%transform sensor data (specifically ax, ay and az) to accel ENU/NED ref
%frame.

%This is to test whether or not it would have a positive impact in the
%overall response of the model (kalman filter). 











%% Yielding Yaw angle via quaternion vectors - Test
%Another test.....
%to get yaw angle (the azimuth angle is usually wrt true north and not
%magnetic north), which is rotation about the zaxis wrt quaternions:

%full rotation q = orientation, 
% while q for azimuth/yaw = [cos(psi/2),0,0, sin(psi/2)]

%multiplying the azimuth_quaternion to the complete quaternion orientation
%vector to get the azimuth angle:
%below is a mathematically derived quaternion for azimuth (shown in report
%in section 2.3.2) :denoted as q_a

[a,b,c,d] = parts(orientation);
coeff_quat = [a,b,c,d];

syms psi
q_a = [cos(psi/2), 0,0, sin(psi/2)]';  %psi is the azimuth/yaw angle need to find
E= coeff_quat(1632,:)*q_a;
s=solve(E==0,psi);
setsolpsi = vpa(s);  %converts solution(s) to numeric value 
degsetsolpsi = vpa(rad2deg(setsolpsi));    %converts rad value to numeric degree

%since two solutions are presented every time from the rad2deg class due to
%the nature of trig functions, will make sure to only select the values
%within the boundary condition of: Real values between -180 deg and +180
%deg. [-180....0....+180]. This is so that orientation and angle are
%preserved. It will be easier to understand which quadrant and at which
%angle the phone's forward axis, in this case Y, is pointing. It'll tell me
%whether or not i'm moving cw or ccw wrt the MN as well. 
%need to do more testing...
R = double(degsetsolpsi);
fpsi = R(R>=-180 & R<=180)


    % for i = 1:length(coeff_quat)
    %     syms psi
    %     q_a = [cos(psi/2), 0,0, sin(psi/2)]';  %psi is the azimuth/yaw angle need to find
    %     E(i,:)= coeff_quat(i,:)*q_a;
    %     s(i,:)=solve(E(i,:)==0,psi);
    %     setsolpsi(i,:) = vpa(s(i,:));  %converts solution(s) to numeric value 
    %     degsetsolpsi(i,:) = vpa(rad2deg(setsolpsi(i,:)));    %converts rad value to numeric degree
    %[s.psi]

    %     R(i,:) = double(degsetsolpsi(i,:));
    %     fpsi(i,:) = R(R(i,:)>=-180 & R(i,:)<=180)
%end 






% if -180 <= degsetsolpsi <=180
%     solpsi = degsetsolpsi;
%     fprintf('%4.2f\n', solpsi)
% end

%At datapoint 1632 --> 109.5 deg
%At datapoint 2715 --> -73.6 deg
%At datapoint 4399 --> -160.11 deg
%At datapoint 5633 --> 42.76 -0.00i deg

%Not dependable approach, however, the values are quite close to the
%eulerAnglesDeg values; in which the rotational offset has been fixed.
%Some of the values may be off due to not taking into account the correct
%directional cosine terms, within the q_a vector, for each
%quaternion based orientation.

%Thus, will end up using the eulerAnglesDeg values instead. 


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

%% Grabbing Azimuth via Azimuth Class - Test

%finding the azimuth, in angles, between two points (Testing phase)
lla62 = lla(62,:);   %capturing the 62nd row of lla table
lla277 = lla(277,:);   % ''
lla384 = lla(384,:);   %capturing the 384th row of lla table
lla516= lla(516,:);    % ''

lat0 = lla0(1,1);
lon0 = lla0(1,2);

lat1 = lla62(1,1);
lon1 = lla62(1,2);

lat2 = lla277(1,1);
lon2 = lla277(1,2);

lat3 = lla384(1,1);
lon3 = lla384(1,2);

lat4 = lla516(1,1);
lon4 = lla516(1,2);

az = azimuth('rh',lat0,lon0,lat3,lon3,'degrees')

%ok, so this is actually placing the starting point at the center of the
%cardinal compass,and the second point is going to be in relation to it.
%Where 0 deg = N, 90 deg = E, 180 deg = South, 270 deg = W and 360 deg = N.

%Is not useful esp, since we want to know the yaw angle wrt orientation at 
%every point, including the starting point. Basically need to know which
%direction the forward axis (in this case,Y) is pointing away from MN
%(magnetic north).

%Thus, lets try using the atan2 class, so as to get rid of determining the
%angle with respect to an initial manually assigned point, as well as, to
%determine where the y-axis is pointing wrt to magnetic north, instead of
%where the device is located wrt the quadrant it's in. 

%% Grabbing Yaw angles from atan2 of the Magnetometer Data - Test

%testing one value at a time...
mag_y = Mag_y(:,1);
mag_x = Mag_x(:,1);

tanaz = atan2(mag_y,mag_x);
tanazd = rad2deg(tanaz);

%datapoint 1632 --> -90 deg
%At datapoint 2715 --> 13.66 deg
%At datapoint 4399 --> 83.86 deg
%At datapoint 5633 --> 174 deg

%This method was treating the E and W as the 0 deg and +/-180 deg axis, 
% respectively. N =+90 deg, S = -90 deg

%According to atan2 rules, the reference axis, or the
%0 deg, is the positive x axis. cw angles are negative, while ccw angles
%are positive. It prefers to travel the shortest route from the positive
%x axis to pi or -pi (or negative x axis). 

%However, atan2 is preferred over atan b/c the latter only showcases angles
%in the 1st and fourth quadrants; where x axis is 0 deg, and positive y is
%90 deg...etc. atan2, however, gives angular values in all quadrants. 
%Yet, despite that, it'll give angles wrt the shortest path from the
%reference axis; which in this case is positive x axis. 

%thus, if x and y are positive values =  quad I. (positive angle)
% if x is negative and y is positive = quad II. (positive angle)
%if x and y are negative values = quad III. (negative angle)
%if x is positive and y is negative = quad IV. (negative angle)

%issue with this is that it would give a quad I angle value when the
%phone's y-axis was pointing W. Ideally it should give a quad III or 
%quad IV read. It reads those values b/c Y and X values were positive (ex:
%datapoint 2715).

%Another issue with this method (of using atan/atan2) is the point of 
%singularity where x = 0. But that is taken care of, mathematically of
%course. 

%it yields angles wrt to the positive or negative values are of x and y, 
%not wrt where y-axis is pointing wrt magnetic north. 
%% Uploading Orientation to get 'Azimuth' angle wrt it - Test

%Load and read Orientation in motion data
Orientation = readtable('Orientation_mobile.xlsx');

%Listing the data based on values and headers
AllOrient = [Orientation.X, Orientation.Y, Orientation.Z];

%Azimuth readings (it's labeled as such in the sensor log of the phone)
OrientZ = AllOrient(:,1);

%All values negative in the list OrientX will be added by 180 deg, while
% positive values will be subtracted by 180 deg. This way, cw values, from
%Magnetic north of 0 deg to south that is 180 deg. Then, values become less
%negative as it approaches MN's 0 deg. This way, + angles denote cw from
%MN to MS while - angles is ccw from 0 deg MN to MS. 

% for i=1:length(OrientX)
%     if OrientX(i,1) <= 0
%         OrientX1(i,1) = abs(OrientX(i,1));
%       
%         
%     elseif OrientX(i,1) >= 0
%         OrientX1(i,1) = OrientX(i,1) - 180;    
% 
%     end  
%   
% end

%This idea was scrapped b/c the Orientation sensor depended on only where
% the forward axis (of the device) was pointing towards. It didin't care
% about magnetic poles. It is very undependable. 
%the orientation sensor is a fusion of the accelerometer and magnetometer.


%% Finding qtrue & qfix;comparing with AHRSFILTER: Yielding Yaw wrt True North - Test
%Let's try another test presented by Matlab devs:
matfile = 'sensorlog_everythang.mat';
SampleRate = 10; % This must match the data rate (Hz) of the phone.

%Calling the ProcessPhoneData function
[Accel, Gyro, Magnet, EulAng] = ProcessPhoneData(matfile);

%based on the function, they're maintaining the order of reading and
%printing data, which is : X-Y-Z, regardless of the reference frame of
%device.
%Thus, Yaw angle would be the last column of qtrue
%printing in X-Y-Z wrt the function's (ProcessPhoneData) eulAngs order (ZXY)

%However, need to adjust the EulAng via interpolation....thus,

qTrue = quaternion([EulAng(:,2), -EulAng(:,3), EulAng(:,1)],'eulerd', 'ZXY', 'frame');

% Get a starting guess at orientation using ecompass. No coefficients
% required. Use the initial orientation estimates to figure out what the
% phone's rotational offset is.
%ecompass requires equal sizes of its argument

Accel_adjusted = AllAccelreads_inter;
Magnet_adjusted = Allmagreads_inter;
q = ecompass(Accel_adjusted, Magnet_adjusted);


Navg = 4; %element size of quaternion
qfix = meanrot(q(1:Navg))./meanrot(qTrue(1:Navg));

% Rotationally corrected phone data via incorporating Declination.
Orientation3 = qfix*qTrue;
eulerAnglesDeg = eulerd(Orientation3,'ZXY','frame');

%This method treats N = 0deg, E = 90 deg, S=+/-180 deg, W = -90deg

%The angular outputs between eulerAnglesDegrees and eulerAnglesDeg have
%about 20 deg difference. The latter has a rotational fix. 

%datapoint 1632 --> 147 deg
%At datapoint 2715 --> -65 deg
%At datapoint 4399 --> -30 deg
%At datapoint 5633 --> 29 deg

%Most accurate method. Will be using this method. 

%% Mathematical Conversions between quaternions and euler angles - Test
%Mathematical equations from wikipedia
%link: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

%The psi angle using quaternions

%quaternions generated from the ecompass class, but need to grab individual
%coefficients.
%where a = q_0, b = q_1, c = q_2, and d = q_3
[a,b,c,d] = parts(orientation);
coeff_quat = [a,b,c,d];

%atan2(Y,X)
psi_quat = atan2(2.*(a.*d+b.*c),1-2.*(c.^2 + d.^2));
psi_quatd = rad2deg(psi_quat);

%datapoint 1632 --> -9.3 deg
%At datapoint 2715 --> 137.6 deg
%At datapoint 4399 -->  -171.2 deg
%At datapoint 5633 --> -131 deg

%it treats S as 0 deg. Left of S, or QuadII 
% is negative angles. However, when phone's y axis is pointing W, angles
% are positive. 

%N = +/-180 deg. Where cw is negative and S = 0 deg. E = -90 deg, W=+90 deg
%It doesn't care about quadrant device is at, only where the y-axis is 
%pointing at wrt the poles. 

%Not ideal, for I prefer to have N = 0 deg. Thus, 'ecompass' class is
%ideal. However, these formulae were used to validate the 'ecompass'
%equation. 

%% Calculating the yaw angle wrt quaternion's azimuth's (q_a) elements - Test

%As stated in an earlier section, the azimuth quaternion is defined as q_a.
% where q_a = [cos(psi/2), 0, 0, sin(psi/2)]. 
%The first element is the real part, while the rest are the coefficients
%for i, j, and k respectively. 
%However, the last three terms are ususally multiplied by the directional
%cosine term. 
%Thus, psi would be different if determining via primarily by cos term
%only.
%Thus, will test via both coefficients, seperately.

[a,b,c,d] = parts(orientation);
coeff_quat = [a,b,c,d];

psi_ang = 2*acos(a);
psi_angd = rad2deg(psi_ang);
%No, not dependable...since the angle is always positive and is usually in
%the range of 150 and 180 deg. This could be b/c of the way the hypersphere
%rotates, where it continously rotates and switches between positive and
%negative values...regardless of orientation. The angle in the cosine term
%is also denoted as the total transformation/rotational angle wrt the axes.

psi_ang2 = 2*asin(d);
psi_ang2d = rad2deg(psi_ang2);
%No, not dependable either...it's stating the angle wrt to the planes. 

%Need for both terms to work together to yield the yaw angle. The EulerD
%seems to be most accurate.

%This method would only work if the there's only movement about the z axis.
%However, that is not the case. Thus, the rotational inputs from roll and
%pitch angles effect the overall angle. 

%local body frame (of device)  vs reference frame (fixed ENU)

%% Plotting the Yaw angles of 'ecompass' EulerAngleDeg variable

%13 deg 3 min is the declination wrt where the data was geographically
% collected. It is subtracted due to being + East.
%Thus, 3 min/60 = 05deg.
%Total declination in deg: 13.05 deg
psi_d = eulerAnglesDegrees(:,1) - 13.05;

figure(3);
plot(psi_d)

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
Bx = [0,0;(dt^2).*cos(psi_r), (dt^2).*sin(psi_r)];
C = [0,1];
Q = [9.3958e-05]; %(variance of the process noise wrt 
% estimated states + accelerometer, where accelerometers -x-y have been 
% equated to have equal variances.)
R = 25; %standard deviation/accuracy is +/- 4.9m. Observable/measured 
%state's noise's variance (GPS in either the X or Y axis).
G = [0;dt^2]; %noise coefficient matrix wrt accelerometer's process
%noise.

%To get the Q matrix, need to grab the variance of the static noise data
%from the phone's accelerometer. Thus, as the phone is statically placed
%horizontally on the table, log stream the data for a few minutes, and then
%calculate the variance of that. This was done in the ME_295A_1 script.

%Calculated the R values in ME_295A_1
%Didn't have to b/c it would lead to an overconfident matrix due
%to the GPS's ability 
R_x = [3.136363636126141e-11];
R_y = [3.379113636223254e-09];

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
eig_A = eig(A); %two positive poles, unstable open system

%Since the Kalman filter is not a controller but an observer does it still
%make sense to test for controllability and observability of the system?

%Cont_x = ctrb(A,Bx); couldn't do it due to the size of the B matrix
%The B matrix isn't constant.

%Let's test for observability
Ob = obsv(A,C);
rank_Ob = rank(Ob); %Must equal number of rows to be full rank, 2
%system is fully observable.

%Since the eigenvalues, E, from the dlqe function are unstable, will use
%the acker observer method.

 

%% Running the Discontinuous Observer Model/Kalman Filter
%Using the interpolated accelerometer values as inputs for filter
a_xx = Accelerometer_inter_x;
a_yy = Accelerometer_inter_y;

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
    x_hat(:,k+1) = A*x_hat(:,k) + B_x*u+L*(z-C*x_hat(:,k));
    
end

figure(5); plot([1:k_end],x_hat)
xlabel('Time step'); ylabel('States')

figure(6); plot([1:k_end],X_GPS)
xlabel('Time step'); ylabel('States')

figure(7); plot([1:k_end], x_hat, [1:k_end], X_GPS)
xlabel('Time step'); ylabel('Comparison btwn Estim. and Meas. States')


%figure(6)
%Plotting the estimated and measured states wrt X
%plot(t,xhat, t,X_GPS)

%Plotting the estimated and measured states wrt Y
%plot(t,yhat, t, Y_GPS)

%% Converting necessary workspace variables to external excel files - Test
%alternative to this process is to save workspace.

%Accelerometer inputs for kalman filter
filename0 = 'a_xx.xls';
writematrix(a_xx,filename0)

filename = 'a_yy.xls';
writematrix(a_yy,filename)


%X_GPS input for kalman filter
filename1 = 'X_GPS.xls';
writematrix(X_GPS,filename1)

%psi_r input for kalman filter
filename2 = 'psi_r.xls';
writematrix(psi_r,filename2)
