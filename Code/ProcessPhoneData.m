function [accel, gyro, mag, eulAngs] = ProcessPhoneData(matfile)
%EXAMPLEHELPERPROCESSPHONEDATA Read logged phone sensor data and convert
%   to the NED frame.
%   Copyright 2020 The MathWorks, Inc.
ld = load(matfile);
% Read Accelerometer data.
accel = [ld.Acceleration.X ld.Acceleration.Y ld.Acceleration.Z];
    
% Read Gyroscope data.
gyro = [ld.AngularVelocity.X ld.AngularVelocity.Y ld.AngularVelocity.Z];
% Read Magnetometer data.
mag =  [ld.MagneticField.X ld.MagneticField.Y ld.MagneticField.Z];
% Read orientation data. [yaw, pitch, roll] based on device ref frame
eulAngs = [ld.Orientation.Z, ld.Orientation.X, ld.Orientation.Y];
% Make sure there are the same number of elements in each sensor array.
na = size(accel,1);
ng = size(gyro,1);
nm = size(mag,1);
