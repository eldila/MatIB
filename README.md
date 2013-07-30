The immersed boundary (IB) method is a mathematical framework for studying fluid-structure interaction initially developed by Charles Peskin to study blood flow through a heart valve [2]. Since its conception, the IB method has found a wide variety of applications in biofluid mechanics and has evolved into a generalized framework [3] for studying fluid-structure interaction problems.

MatIB is a simple Matlab implementation of the IB method that allows students and researchers to solve simple fluid-structure interaction problems with minimal overhead. The algorithm employed in MatIB is described in Peskin’s review paper [3] and is an adaptation of the Lai-Peskin algorithm [1]. For clarity, we have limited the scope of our implementation and this code should therefore not be thought as a generalized IB toolkit. Instead, MatIB’s codebase acts as a foundation for further experimentation and extension.

MatIB is released under the MIT Open Source License and is free to use for any purpose. However, any resulting publications should cite the MatIB User Guide.

References:

[1] Ming-Chih Lai and Charles S. Peskin. An immersed boundary method with formal second-order accuracy and reduced numerical viscosity. Journal of Computational Physics, 160(2):705–719, 2000.

[2] Charles S. Peskin. Flow patterns around heart valves: A numerical method. Journal of Computational Physics, 10(2):252 – 271, 1972.

[3] Charles S. Peskin. The immersed boundary method. Acta Numerica, 11:479– 517, 2002.