# Waiter Robot
The following is a MATLAB code to simulate the guidance and navigation of a waiter robot in a japanese restaurant. The project was part of the subject "Guidance and Navigation of Robots" at Technical University of Madrid (UPM). The program makes use of the simulation program Apolo, which communicates with MATLAB. The employed robot is the Pioneer 3-AT with a differential locomotion system. It accounts with a laser telemeter to locate itself in the restaurant thanks to several beacons disposed along the restaurant using Extended Kalman Filter algorithm. The robot additionally includes three ultrasonic sensors to be able to perform react control. The program uses the A* algorithm for free-collision trajectory planning. This repository includes the files used for sensor calibration, localization, planning and control.

[<img src="https://github.com/danisotelo/robot_waiter/blob/main/img/readme_img.png" width="100%">](https://www.youtube.com/watch?v=LhBZU4BIHkA)

## Getting Started
### Cloning the Repository
To clone the repository and start using the **waiter robot**, follow these steps:
1. Open a terminal or command prompt
2. Navigate to the directory where you want to clone the repository.
3. Run the following command:
```
git clone https://github.com/danisotelo/robot_waiter.git
```
### Installing Apolo
Before running the program, you need to install Apolo, which is a simulation environment for Windows developed by ETSIDI (UPM) professors. This software enables the map definition and the graphic visualization of the robot and its sensors in that environment. It also simulates the kinematic model of the robot and the models of the sensor measurements. The program can be downloaded and installed from https://github.com/mhernando/Apolo/tree/master.

Download the Apolo folder in your desktop and add the ```Apolo/Matlab``` folder. Finally, move the MATLAB files you downloaded from this repository into the ```Apolo/Matlab``` folder.

## Running the Program
Open Apolo and load the restaurant environment from the ```restaurant.xml``` file. Next, open MATLAB and run the main program ```navegador_principal.m```. This will start the demonstration of the robot executing the trajectory to carry the food to the desired tables. A video of the demonstrator can be seen at: https://www.youtube.com/watch?v=LhBZU4BIHkA.
