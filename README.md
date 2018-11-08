# Unscented Kalman Filter Project

This repository contains my solution for the project "Unscented Kalman Filter Project" of the Udacity Self-Driving Car Engineer Nanodegree Program. A description of a basic setup can be found in the [original repository](https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project). The written code could be found in the files [ukf.cpp](./src/ukf.cpp), [tools.cpp](./src/tools.cpp) and their corresponding header files.

The following part of the README contains a very short writeup which describes what is done.

---

In this project an Unscented Kalman filter is implemented to combine laser and radar measurements to predict the location of a car.

In the file [ukf.cpp](./src/ukf.cpp) the initialization and logic of the UKF is implemented. On an arriving measurement first the prediction step is preformed and then the update step. If there are both laser and radar measurements at almost the same time the second prediction step is skipped to gain some performance. If the arriving measurement is the first measurement the UKF is initialized instead of preforming a prediction update step.

Also the computation of the prediction and update step is done in this file. These computations simply follow the equations for an UKF.

In [tools.cpp](./src/tools.cpp) there is a function for calculating the error RMSE.

---

Running the implementation with the [simulator](https://github.com/udacity/self-driving-car-sim/releases) leads to the following results for the RMSE:


| RMSE UKF | Dataset 1 Both | Dataset 1 Laser | Dataset 1 Radar | Dataset 2 Both | Dataset 2 Laser | Dataset 2 Radar |
|------|----------------|-----------------|-----------------|----------------|-----------------|-----------------|
| X    | **0.0755**     | 0.1259          | 0.1891          | **0.0959**     | 0.1107          | 0.2032          |
| Y    | **0.0836**     | 0.0982          | 0.2755          | **0.0813**     | 0.0951          | 0.2306          |
| VX   | **0.3162**     | 0.7033          | 0.3834          | **0.4986**     | 0.6530          | 0.5749          |
| VY   | **0.2016**     | 0.2316          | 0.3427          | **0.2267**     | 0.2460          | 0.3460          |

With these results the RMSE of the implementation is low enough to successfully track the car.


As one can see combining laser and radar measurements significantly improve the performance if as one would only use one sensor. Additionally the data shows that the laser measurements provide a way more accurate prediction. While the VX, VY RMSE of the radar sensor is not that bad compared to the laser sensor, the RMSE of X and Y is up to 3 times higher for the radar sensor. The reason for this should be that the position X and Y can only be measured indirectly by the radar sensor while the laser sensor measures it directly.

Recalling the following results from the previous project for an EKF

| RMSE EKF | Dataset 1 Both | Dataset 1 Laser | Dataset 1 Radar | Dataset 2 Both | Dataset 2 Laser | Dataset 2 Radar |
|------|----------------|-----------------|-----------------|----------------|-----------------|-----------------|
| X    | **0.0973**     | 0.1473          | 0.2302          | **0.0726**     | 0.1169          | 0.2706          |
| Y    | **0.0855**     | 0.1153          | 0.3464          | **0.0965**     | 0.1260          | 0.3853          |
| VX   | **0.4513**     | 0.6383          | 0.5835          | **0.4216**     | 0.6227          | 0.6524          |
| VY   | **0.4399**     | 0.5346          | 0.8040          | **0.4932**     | 0.6024          | 0.9218          |

one could see that the UKF is almost every time better then the EKF. The RMSE for X and Y is quite the same for the UKF and EKF but there is a significant improvement for VX and VY.

Calculating the NIS for laser and radar predictions results in

  | x².950  | x².900  | x².100  | x².050
--|---|---|---|--
 Radar NIS | 0.919679  | 0.851406  | 0.0803213  |  0.0401606
 Laser NIS | 0.947791  | 0.911647  | 0.060241  |  0.0281125

which nearly matches the expected distribution. This shows that the process noise is well adjusted.


# Build & Run

### Dependencies

* cmake >= 3.5
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
* gcc/g++ >= 5.4

### Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ExtendedKF `
