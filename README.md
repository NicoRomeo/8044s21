# 8044s21
Simulations and visualizations for 8.044, spring 21 edition

## Hard-disk simulations

We follow Werener Krauth's "Statistical Mechanics and Computations" algorithm to directly slve the exact dynamics of the spheres through an event-driven algorithm: namely, we exactly solve Newton's laws of motion for all particles until the next collision, between disks or betwen a disk and a wall happens. We then interpolate the result to sample the dynamics at regular time points to contrstruct a movie.
