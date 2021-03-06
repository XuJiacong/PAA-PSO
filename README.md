# PID-controller-based Auxiliary Algorithm for PSO

Official implementation of: A PID Controller Architecture Inspired Enhancement to the PSO Algorithm (FICC 2022)

![highlights](figs/analogy.png)

## Highlights
- **PSO Algorithm is an I controller**: The accumulation of the error in PSO algorithm serves the same function as I controller.
- **Overshoot problem of I controller**: The PSO algorithm suffers from long-time convergence due to oscillation around minima.
- **Added D and P controllers to PSO**: To make up a full PID controller with 3 DOF and D component will moderate the overshoot.
- **Improvement**; This method could signicantly accelerate the convergence speed of PSO algorithm and improve the optimization accuracy with the sacrice of time usage of less than 2%.


## Experiments
we applied all the algorithms on the 30 functions in CEC2014 benchmark test suite for real-parameter optimization.

## Accuracy
PBSv2 is another proposed method based on previous work.
![highlights](figs/table.png)

## Convergence
![highlights](figs/func18_all.jpg)

## Reference
Most of our code was developed based on (http://www5.zzu.edu.cn/cilab/Benchmark/wysyhwtcsj.htm), thanks for their nice contribution.
