# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

---

[//]: # (Image References)
[image1]: ./images/plot1.png
[image2]: ./images/plot2.png
[image3]: ./images/plot1_rmse.png
[image4]: ./images/plot2_rmse.png


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Results

[***Path 1***](./data/sample-laser-radar-measurement-data-1.txt)

![alt text][image1]

![alt text][image3]

[***Path 2***](./data/sample-laser-radar-measurement-data-2.txt)

![alt text][image2]

![alt text][image4]
