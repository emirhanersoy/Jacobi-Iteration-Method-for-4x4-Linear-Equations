using System;

public static class JacobiIteration
{
      // EMIRHAN ERSOY
      // 22118080001
    public static (double[], int, double) Solve(double[,] A, double[] b, double tolerance, int maxIterations)
    {
        if (A.GetLength(0) != A.GetLength(1) || A.GetLength(0) != b.Length)
        {
            throw new ArgumentException("Matrix and vector dimensions must be compatible.");
        }

        int n = A.GetLength(0);

        double[] diagonal = new double[n];
        for (int i = 0; i < n; i++)
        {
            diagonal[i] = A[i, i];
        }

        double[] bVector = b; //result vector

        double[,] offDiagonal = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    offDiagonal[i, j] = A[i, j];
                }
            }
        }

        double[] invDiagonal = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (diagonal[i] == 0)
            {
                throw new ArgumentException("Matrix is not diagonally dominant (zero diagonal element).");
            }
            invDiagonal[i] = 1.0 / diagonal[i];
        }

        double[] x = new double[n];
        double[] xOld = new double[n]; // X values for iterations

        double error = double.PositiveInfinity;
        int iterations = 0;

        while (error > tolerance && iterations < maxIterations)
        {
            iterations++;

            Array.Copy(x, xOld, n);

            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                    {
                        sum += offDiagonal[i, j] * xOld[j];
                    }
                }
                x[i] = invDiagonal[i] * (bVector[i] - sum);
            }

            error = 0.0;
            for (int i = 0; i < n; i++)
            {
                error = Math.Max(error, Math.Abs(x[i] - xOld[i]));
            }
        }

        if (error > tolerance)
        {
            Console.WriteLine("Warning: Convergence tolerance was not met at termination.");
        }

        return (x, iterations, error);
    }
    public static void Main(string[] args)
    {
        // Define the matrix A
        double[,] A = new double[,] {
            { 1, 3, 1, -1 },
            { 2, 0, 1, 1},
            { 0, -1, 4, 1 },
            {0, 1, 1, -5 }
        };

        // Define the result vector b
        double[] b = { 1, 3, 6, 16};

        // Set tolerance and maximum iterations
        double tolerance = 1e-6; //generally 1*10^-6 is used
        int maxIterations = 100; 

        // Solve the system using Jacobi iteration
        var (solution, iterations, error) = Solve(A, b, tolerance, maxIterations);

        // Print the solution
        Console.WriteLine("Solution:");
        for (int i = 0; i < solution.Length; i++)
        {
            Console.WriteLine($"x[{i + 1}] = {solution[i]}");
        }

        Console.WriteLine($"\nIterations: {iterations}");
        Console.WriteLine($"Error: {error}");
        Console.ReadKey();
    }
}
