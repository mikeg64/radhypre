minimal working example of solving a linear system Ax=b using HYPRE’s IJ interface. 

This will demonstrate the basics: setting up the matrix, right-hand side, solution vector, and calling a solver (BoomerAMG in this case). 




Solution: x = [0.090909, 0.636364]

which is the solution to:
\left[ \begin{matrix}4&1\\ 1&3\end{matrix}\right] \left[ \begin{matrix}x_1\\ x_2\end{matrix}\right] =\left[ \begin{matrix}1\\ 2\end{matrix}\right] 



mkdir build && cd build
cmake ..
make
mpirun -np 1 ./hypre_example




