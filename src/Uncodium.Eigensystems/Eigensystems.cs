using Aardvark.Base;
using System;
using System.Runtime.CompilerServices;

namespace Uncodium
{
    /// <summary>
    /// Numerical diagonalization of 3x3 matrices.
    /// See original work at https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/.
    /// 
    /// A common scientific problem is the numerical calculation of the
    /// eigensystem of symmetric or hermitian 3x3 matrices. If this
    /// calculation has to be performed many times, standard packages
    /// like LAPACK, the GNU Scientific Library, and the Numerical Recipes
    /// Library may not be the optimal choice because they are optimized
    /// mainly for large matrices.
    /// 
    /// This is a C# port of the C implementations of several algorithms
    /// which were specifically optimized for 3x3 problems.
    /// All algorithms are discussed in detail in the following paper:
    /// 
    /// Joachim Kopp
    /// Efficient numerical diagonalization of hermitian 3x3 matrices
    /// Int.J.Mod.Phys.C 19 (2008) 523-548
    /// arXiv.org: http://arxiv.org/abs/physics/0610206
    /// </summary>
    public static class Eigensystems
    {
        /// <summary>
        /// Calculates the eigensystem of a real symmetric 2x2 matrix
        ///   [ A  B ]
        ///   [ B  C ]
        /// in the form
        ///   [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
        ///   [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
        /// where rt1 >= rt2. Note that this convention is different from the one used
        /// in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
        /// </summary>
        public static void Dsyev2(this in M22d A, out M22d Q, out V2d w)
        {
            Dsyev2(A.M00, A.M01, A.M11, out double rt1, out double rt2, out double cs, out double sn);
            Q = new M22d(cs, sn, -sn, cs);
            w = new V2d(rt1, rt2);
        }

        /// <summary>
        /// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
        /// analytical algorithm. Only the diagonal and upper triangular parts of A
        /// are accessed. The access is read-only.
        /// </summary>
        /// <param name="A">The symmetric input matrix.</param>
        /// <param name="w">Eigenvalues.</param>
        /// <returns>Eigenvalues of matrix A.</returns>
        public static void Dsyevc3(this in M33d A, out V3d w)
        {
            // Determine coefficients of characteristic poynomial. We write
            //       | a   d   f  |
            //  A =  | d*  b   e  |
            //       | f*  e*  c  |
            double a = A.M00, b = A.M11, c = A.M22, d = A.M01, e = A.M12, f = A.M02;
            double de = d * e;  // d * e
            double dd = d * d;  // d^2
            double ee = e * e;  // e^2
            double ff = f * f;  // f^2
            var m = a + b + c;
            // a*b + a*c + b*c - d^2 - e^2 - f^2
            var c1 = (a * b + a * c + b * c) - (dd + ee + ff);
            // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)
            var c0 = c * dd + a * ee + b * ff - a * b * c - 2.0 * f * de;

            var p = m * m - 3.0 * c1;
            var q = m * (p - (3.0 / 2.0) * c1) - (27.0 / 2.0) * c0;
            var sqrt_p = Math.Sqrt(Math.Abs(p));

            var phi = 27.0 * (0.25 * c1 * c1 * (p - c1) + c0 * (q + 27.0 / 4.0 * c0));
            phi = (1.0 / 3.0) * Math.Atan2(Math.Sqrt(Math.Abs(phi)), q);

            var cos = sqrt_p * Math.Cos(phi);
            var sin = (1.0 / M_SQRT3) * sqrt_p * Math.Sin(phi);

            var wy = (1.0 / 3.0) * (m - cos);
            var wz = wy + sin;
            var wx = wy + cos;
            wy -= sin;

            w = new V3d(wx, wy, wz);
        }
        
        /// <summary>
        /// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
        /// matrix A using Cuppen's Divide and Conquer algorithm.
        /// The function accesses only the diagonal and upper triangular parts of A.
        /// The access is read-only.
        /// </summary>
        /// <param name="A">The symmetric input matrix.</param>
        /// <param name="Q">Eigenvectors.</param>
        /// <param name="w">Eigenvalues.</param>
        public static void Dsyevd3(this in M33d A, out M33d Q, out V3d w)
        {
            //Dsyevd3(Unpack(A), out double[,] _Q, out double[] _w);
            //Q = Pack(_Q);
            //w = Pack(_w);
            M33d R;                // Householder transformation matrix
            var P = M33d.Zero;   // Unitary transformation matrix which diagonalizes D + w w^T
            V2d e;      // Off-diagonal elements after Householder transformation
            V3d d;      // Eigenvalues of split matrix in the "divide" step)
            double c, s;                // Eigenvector of 2x2 block in the "divide" step
            V3d z;      // Numerators of secular equation / Updating vector
            double t;                   // Miscellaenous temporary stuff

            // Transform A to real tridiagonal form by the Householder method.
            Dsytrd3(A, out R, out w, out e);

            // "Divide"
            // --------

            // Detect matrices that factorize to avoid multiple eigenvalues in the Divide/Conquer algorithm.
            for (int i = 0; i < 2; i++)
            {
                t = Math.Abs(w[i]) + Math.Abs(w[i + 1]);
                if (Math.Abs(e[i]) <= 8.0 * DBL_EPSILON * t)
                {
                    if (i == 0)
                    {
                        Dsyev2(w.Y, e.Y, w.Z, out d.Y, out d.Z, out c, out s);
                        w.Y = d.Y;
                        w.Z = d.Z;

                        Q = new M33d(
                            1.0, 0.0, 0.0,
                            0.0, s * R.M12 + c * R.M11, c * R.M12 - s * R.M11,
                            0.0, s * R.M22 + c * R.M21, c * R.M22 - s * R.M21
                            );
                    }
                    else
                    {
                        Dsyev2(w[0], e[0], w[1], out d.X, out d.Y, out c, out s);
                        w.X = d.X;
                        w.Y = d.Y;

                        Q = new M33d(
                            c, -s, 0,
                            R.M11 * s, R.M11 * c, R.M12,
                            R.M21 * s, R.M21 * c, R.M22
                            );
                    }

                    return;
                }
            }

            // Calculate eigenvalues and eigenvectors of 2x2 block.
            Dsyev2(w.Y - e.X, e.Y, w.Z, out d.Y, out d.Z, out c, out s);
            d.X = w.X - e.X;

            // "Conquer"
            // ---------

            // Determine coefficients of secular equation.
            z.X = e.X;
            z.Y = e.X * c * c;
            z.Z = e.X * s * s;

            // Call slvsec3 with d sorted in ascending order. We make
            // use of the fact that dsyev2 guarantees d[1] >= d[2].
            if (d.X < d.Z)
                Slvsec3(d, z, out w, out P, 0, 2, 1);
            else if (d.X < d.Y)
                Slvsec3(d, z, out w, out P, 2, 0, 1);
            else
                Slvsec3(d, z, out w, out P, 2, 1, 0);

            // Calculate eigenvectors of matrix D + beta * z * z^t and store them in the columns of P.
            z.X = Math.Sqrt(Math.Abs(e.X));
            z.Y = c * z.X;
            z.Z = -s * z.X;

            // Detect duplicate elements in d to avoid division by zero.
            t = 8.0 * DBL_EPSILON * (Math.Abs(d.X) + Math.Abs(d.Y) + Math.Abs(d.Z));
            if (Math.Abs(d.Y - d.X) <= t)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (P[0, j] * P[1, j] <= 0.0)
                    {
                        P[0, j] = z.Y;
                        P[1, j] = -z.X;
                        P[2, j] = 0.0;
                    }
                    else
                    {
                        for (int i = 0; i < 3; i++)
                            P[i, j] = z[i] / P[i, j];
                    }
                }
            }
            else if (Math.Abs(d.Z - d.X) <= t)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (P[0, j] * P[2, j] <= 0.0)
                    {
                        P[0, j] = z.Z;
                        P[1, j] = 0.0;
                        P[2, j] = -z.X;
                    }
                    else
                    {
                        for (int i = 0; i < 3; i++)
                            P[i, j] = z[i] / P[i, j];
                    }
                }
            }
            else
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        if (P[i, j] == 0.0)
                        {
                            P[i, j] = 1.0;
                            P[(i + 1) % 3, j] = 0.0;
                            P[(i + 2) % 3, j] = 0.0;
                            break;
                        }
                        else
                        {
                            P[i, j] = z[i] / P[i, j];
                        }
                    }
                }
            }

            // Normalize eigenvectors of D + beta * z * z^t.
            for (int j = 0; j < 3; j++)
            {
                t = P[0, j] * P[0, j] + P[1, j] * P[1, j] + P[2, j] * P[2, j];
                t = 1.0 / Math.Sqrt(t);
                for (int i = 0; i < 3; i++)
                    P[i, j] *= t;
            }

            // Undo diagonalization of 2x2 block.
            for (int j = 0; j < 3; j++)
            {
                t = P[1, j];
                P[1, j] = c * t - s * P[2, j];
                P[2, j] = s * t + c * P[2, j];
            }

            // Undo Householder transformation.
            Q = M33d.Zero;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    t = P[k, j];
                    for (int i = 0; i < 3; i++)
                    {
                        Q[i, j] += t * R[i, k];
                    }
                }
            }
        }

        /// <summary>
        /// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
        /// matrix A using Cardano's method for the eigenvalues and an analytical
        /// method based on vector cross products for the eigenvectors. However,
        /// if conditions are such that a large error in the results is to be
        /// expected, the routine falls back to using the slower, but more
        /// accurate QL algorithm. Only the diagonal and upper triangular parts of A need
        /// to contain meaningful values. Access to A is read-only.
        /// </summary>
        /// <param name="A">The symmetric input matrix.</param>
        /// <param name="Q">Eigenvectors.</param>
        /// <param name="w">Eigenvalues.</param>
        /// <returns>True (success), false (no convergence).</returns>
        public static bool Dsyevh3(this in M33d A, out M33d Q, out V3d w)
        {
            double norm;          // Squared norm or inverse norm of current eigenvector
                                  //  double n0, n1;        // Norm of first and second columns of A
            double error;         // Estimated maximum roundoff error
            double t, u;          // Intermediate storage

            // Calculate eigenvalues
            Dsyevc3(A, out w);

            //  n0 = SQR(A[0][0]) + SQR(A[0][1]) + SQR(A[0][2]);
            //  n1 = SQR(A[0][1]) + SQR(A[1][1]) + SQR(A[1][2]);

            t = Math.Abs(w.X);
            if ((u = Math.Abs(w.Y)) > t)
                t = u;
            if ((u = Math.Abs(w.Z)) > t)
                t = u;
            if (t < 1.0)
                u = t;
            else
                u = t * t;
            error = 256.0 * DBL_EPSILON * u * u;
            //  error = 256.0 * DBL_EPSILON * (n0 + u) * (n1 + u);

            Q = new M33d(
                0.0, A.M01 * A.M12 - A.M02 * A.M11, 0.0,
                0.0, A.M02 * A.M01 - A.M12 * A.M00, 0.0,
                0.0, A.M01 * A.M01, 0.0
                );

            // Calculate first eigenvector by the formula
            //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
            Q.M00 = Q.M01 + A.M02 * w.X;
            Q.M10 = Q.M11 + A.M12 * w.X;
            Q.M20 = (A.M00 - w.X) * (A.M11 - w.X) - Q.M21;
            norm = Q.M00 * Q.M00 + Q.M10 * Q.M10 + Q.M20 * Q.M20;

            // If vectors are nearly linearly dependent, or if there might have
            // been large cancellations in the calculation of A[i][i] - w[0], fall
            // back to QL algorithm
            // Note that this simultaneously ensures that multiple eigenvalues do
            // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
            // i.e. all columns of A - w[0] * I are linearly dependent.
            if (norm <= error)
                return Dsyevq3(A, out Q, out w);
            else                      // This is the standard branch
            {
                norm = Math.Sqrt(1.0 / norm);
                Q.M00 = Q.M00 * norm;
                Q.M10 = Q.M10 * norm;
                Q.M20 = Q.M20 * norm;

            }

            // Calculate second eigenvector by the formula
            //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
            Q.M01 = Q.M01 + A.M02 * w.Y;
            Q.M11 = Q.M11 + A.M12 * w.Y;
            Q.M21 = (A.M00 - w.Y) * (A.M11 - w.Y) - Q.M21;
            norm = Q.M01 * Q.M01 + Q.M11 * Q.M11 + Q.M21 * Q.M21;
            if (norm <= error)
                return Dsyevq3(A, out Q, out w);
            else
            {
                norm = Math.Sqrt(1.0 / norm);
                Q.M01 = Q.M01 * norm;
                Q.M11 = Q.M11 * norm;
                Q.M21 = Q.M21 * norm;
            }

            // Calculate third eigenvector according to
            //   v[2] = v[0] x v[1]
            Q.M02 = Q.M10 * Q.M21 - Q.M20 * Q.M11;
            Q.M12 = Q.M20 * Q.M01 - Q.M00 * Q.M21;
            Q.M22 = Q.M00 * Q.M11 - Q.M10 * Q.M01;

            return true;
        }
        /// <summary>
        /// Same as Dsyevh3, but additional order of Eigenvalues (ascending).
        /// </summary>
        public static bool Dsyevh3asc(this in M33d A, out M33d Q, out V3d w, out int[] order)
        {
            if (Dsyevh3(A, out Q, out w))
            {
                order = GetOrderAscending(in w);
                return true;
            }
            else
            {
                order = s_order012;
                return false;
            }
        }
        /// <summary>
        /// Same as Dsyevh3, but additional order of Eigenvalues (descending).
        /// </summary>
        public static bool Dsyevh3desc(this in M33d A, out M33d Q, out V3d w, out int[] order)
        {
            if (Dsyevh3(A, out Q, out w))
            {
                order = GetOrderDescending(in w);
                return true;
            }
            else
            {
                order = s_order012;
                return false;
            }
        }

        private static int[] GetOrderAscending(in V3d w)
        {
            if (w.X < w.Y)
            {
                if (w.X < w.Z)
                {
                    return w.Y < w.Z ? s_order012 : s_order021;
                }
                else
                {
                    return s_order201;
                }
            }
            else
            {
                if (w.X < w.Z)
                {
                    return s_order102;
                }
                else
                {
                    return w.Y < w.Z ? s_order120 : s_order210;
                }
            }
        }
        private static int[] GetOrderDescending(in V3d w)
        {
            if (w.X > w.Y)
            {
                if (w.X > w.Z)
                {
                    return w.Y > w.Z ? s_order012 : s_order021;
                }
                else
                {
                    return s_order201;
                }
            }
            else
            {
                if (w.X > w.Z)
                {
                    return s_order102;
                }
                else
                {
                    return w.Y > w.Z ? s_order120 : s_order210;
                }
            }
        }
        private static readonly int[] s_order012 = new[] { 0, 1, 2 };
        private static readonly int[] s_order021 = new[] { 0, 2, 1 };
        private static readonly int[] s_order102 = new[] { 1, 0, 2 };
        private static readonly int[] s_order120 = new[] { 1, 2, 0 };
        private static readonly int[] s_order201 = new[] { 2, 0, 1 };
        private static readonly int[] s_order210 = new[] { 2, 1, 0 };

        /// <summary>
        /// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
        /// matrix A using the Jacobi algorithm.
        /// The upper triangular part of A is destroyed during the calculation,
        /// the diagonal elements are read but not destroyed, and the lower
        /// triangular elements are not referenced at all.
        /// </summary>
        /// <param name="A">The symmetric input matrix.</param>
        /// <param name="Q">Eigenvectors.</param>
        /// <param name="w">Eigenvalues.</param>
        /// <returns>True (success), false (no convergence).</returns>
        public static bool Dsyevj3(this ref M33d A, out M33d Q, out V3d w)
        {
            double sd, so;                  // Sums of diagonal resp. off-diagonal elements
            double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
            double g, h, z, theta;          // More temporary storage
            double thresh;

            // Initialize Q to the identitity matrix
            Q = M33d.Identity;

            // Initialize w to diag(A)
            w = new V3d(A.M00, A.M11, A.M22);

            // Calculate SQR(tr(A))  
            sd = 0.0;
            sd += Math.Abs(w.X);
            sd += Math.Abs(w.Y);
            sd += Math.Abs(w.Z);
            sd = sd * sd;

            // Main iteration loop
            for (int nIter = 0; nIter < 50; nIter++)
            {
                // Test for convergence 
                so = 0.0;
                for (int p = 0; p < 3; p++)
                    for (int q = p + 1; q < 3; q++)
                        so += Math.Abs(A[p, q]);
                if (so == 0.0)
                    return true;

                if (nIter < 4)
                    thresh = 0.2 * so / (3 * 3);
                else
                    thresh = 0.0;

                // Do sweep
                for (int p = 0; p < 3; p++)
                    for (int q = p + 1; q < 3; q++)
                    {
                        g = 100.0 * Math.Abs(A[p, q]);
                        if (nIter > 4 && Math.Abs(w[p]) + g == Math.Abs(w[p])
                                       && Math.Abs(w[q]) + g == Math.Abs(w[q]))
                        {
                            A[p, q] = 0.0;
                        }
                        else if (Math.Abs(A[p, q]) > thresh)
                        {
                            // Calculate Jacobi transformation
                            h = w[q] - w[p];
                            if (Math.Abs(h) + g == Math.Abs(h))
                            {
                                t = A[p, q] / h;
                            }
                            else
                            {
                                theta = 0.5 * h / A[p, q];
                                if (theta < 0.0)
                                    t = -1.0 / (Math.Sqrt(1.0 + theta * theta) - theta);
                                else
                                    t = 1.0 / (Math.Sqrt(1.0 + theta * theta) + theta);
                            }
                            c = 1.0 / Math.Sqrt(1.0 + t * t);
                            s = t * c;
                            z = t * A[p, q];

                            // Apply Jacobi transformation
                            A[p, q] = 0.0;
                            w[p] -= z;
                            w[q] += z;
                            for (int r = 0; r < p; r++)
                            {
                                t = A[r, p];
                                A[r, p] = c * t - s * A[r, q];
                                A[r, q] = s * t + c * A[r, q];
                            }
                            for (int r = p + 1; r < q; r++)
                            {
                                t = A[p, r];
                                A[p, r] = c * t - s * A[r, q];
                                A[r, q] = s * t + c * A[r, q];
                            }
                            for (int r = q + 1; r < 3; r++)
                            {
                                t = A[p, r];
                                A[p, r] = c * t - s * A[q, r];
                                A[q, r] = s * t + c * A[q, r];
                            }

                            // Update eigenvectors        
                            for (int r = 0; r < 3; r++)
                            {
                                t = Q[r, p];
                                Q[r, p] = c * t - s * Q[r, q];
                                Q[r, q] = s * t + c * Q[r, q];
                            }
                        }
                    }
            }

            return false;
        }

        /// <summary>
        /// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
        /// matrix A using the QL algorithm with implicit shifts, preceded by a
        /// Householder reduction to tridiagonal form.
        /// The function accesses only the diagonal and upper triangular parts of A.
        /// The access is read-only.
        /// </summary>
        /// <param name="A">The symmetric input matrix.</param>
        /// <param name="Q">Eigenvectors.</param>
        /// <param name="w">Eigenvalues.</param>
        /// <returns>True (success), false (no convergence).</returns>
        public static bool Dsyevq3(this in M33d A, out M33d Q, out V3d w)
        {
            var e = V3d.Zero;               // The third element is used only as temporary workspace
            double g, r, p, f, b, s, c, t;  // Intermediate storage
            int nIter;
            int m;

            // Transform A to real tridiagonal form by the Householder method
            Dsytrd3(A, out Q, out w, out V2d _e);
            e.X = _e.X; e.Y = _e.Y;

            // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
            // with the QL method
            //
            // Loop over all off-diagonal elements
            for (int l = 0; l < 2; l++)
            {
                nIter = 0;
                while (true)
                {
                    // Check for convergence and exit iteration loop if off-diagonal
                    // element e(l) is zero
                    for (m = l; m <= 1; m++)
                    {
                        g = Math.Abs(w[m]) + Math.Abs(w[m + 1]);
                        if (Math.Abs(e[m]) + g == g)
                            break;
                    }
                    if (m == l)
                        break;

                    if (nIter++ >= 30)
                        return false;

                    // Calculate g = d_m - k
                    g = (w[l + 1] - w[l]) / (e[l] + e[l]);
                    r = Math.Sqrt(g * g + 1.0);
                    if (g > 0)
                        g = w[m] - w[l] + e[l] / (g + r);
                    else
                        g = w[m] - w[l] + e[l] / (g - r);

                    s = c = 1.0;
                    p = 0.0;
                    for (int i = m - 1; i >= l; i--)
                    {
                        f = s * e[i];
                        b = c * e[i];
                        if (Math.Abs(f) > Math.Abs(g))
                        {
                            c = g / f;
                            r = Math.Sqrt(c * c + 1.0);
                            e[i + 1] = f * r;
                            c *= (s = 1.0 / r);
                        }
                        else
                        {
                            s = f / g;
                            r = Math.Sqrt(s * s + 1.0);
                            e[i + 1] = g * r;
                            s *= (c = 1.0 / r);
                        }

                        g = w[i + 1] - p;
                        r = (w[i] - g) * s + 2.0 * c * b;
                        p = s * r;
                        w[i + 1] = g + p;
                        g = c * r - b;

                        // Form eigenvectors
                        for (int k = 0; k < 3; k++)
                        {
                            t = Q[k, i + 1];
                            Q[k, i + 1] = s * Q[k, i] + c * t;
                            Q[k, i] = c * Q[k, i] - s * t;
                        }
                    }
                    w[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            }

            return true;
        }

        /// <summary>
        /// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
        /// matrix A using Cardano's method for the eigenvalues and an analytical
        /// method based on vector cross products for the eigenvectors.
        /// Only the diagonal and upper triangular parts of A need to contain meaningful
        /// values. However, all of A may be used as temporary storage and may hence be
        /// destroyed.
        /// </summary>
        /// <param name="A">The symmetric 3x3 input matrix.</param>
        /// <param name="Q">Eigenvectors.</param>
        /// <param name="w">Eigenvalues.</param>
        public static void Dsyevv3(this ref M33d A, out M33d Q, out V3d w)
        {
            double norm;          // Squared norm or inverse norm of current eigenvector
            double n0, n1;        // Norm of first and second columns of A
            double n0tmp, n1tmp;  // "Templates" for the calculation of n0/n1 - saves a few FLOPS
            double thresh;        // Small number used as threshold for floating point comparisons
            double error;         // Estimated maximum roundoff error in some steps
            double wmax;          // The eigenvalue of maximum modulus
            double f, t;          // Intermediate storage
            int i, j;             // Loop counters
            Q = M33d.Zero;

            // Calculate eigenvalues
            Dsyevc3(A, out w);

            wmax = Math.Abs(w.X);
            if ((t = Math.Abs(w.Y)) > wmax)
                wmax = t;
            if ((t = Math.Abs(w.Z)) > wmax)
                wmax = t;
            thresh = (8.0 * DBL_EPSILON * wmax) * (8.0 * DBL_EPSILON * wmax);

            // Prepare calculation of eigenvectors
            n0tmp = A.M01 * A.M01 + A.M02 * A.M02;
            n1tmp = A.M01 * A.M01 + A.M12 * A.M12;
            Q.M01 = A.M01 * A.M12 - A.M02 * A.M11;
            Q.M11 = A.M02 * A.M01 - A.M12 * A.M00;
            Q.M21 = A.M01 * A.M01;

            // Calculate first eigenvector by the formula
            //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
            A.M00 -= w.X;
            A.M11 -= w.X;
            Q.M00 = Q.M01 + A.M02 * w.X;
            Q.M10 = Q.M11 + A.M12 * w.X;
            Q.M20 = A.M00 * A.M11 - Q.M21;
            norm = Q.M00 * Q.M00 + Q.M10 * Q.M10 + Q.M20 * Q.M20;
            n0 = n0tmp + A.M00 * A.M00;
            n1 = n1tmp + A.M11 * A.M11;
            error = n0 * n1;

            if (n0 <= thresh)         // If the first column is zero, then (1,0,0) is an eigenvector
            {
                Q.M00 = 1.0;
                Q.M10 = 0.0;
                Q.M20 = 0.0;
            }
            else if (n1 <= thresh)    // If the second column is zero, then (0,1,0) is an eigenvector
            {
                Q.M00 = 0.0;
                Q.M10 = 1.0;
                Q.M20 = 0.0;
            }
            else if (norm < SQR(64.0 * DBL_EPSILON) * error)
            {                         // If angle between A[0] and A[1] is too small, don't use
                t = SQR(A.M01);       // cross product, but calculate v ~ (1, -A0/A1, 0)
                f = -A.M00 / A.M01;
                if (SQR(A.M11) > t)
                {
                    t = SQR(A.M11);
                    f = -A.M01 / A.M11;
                }
                if (SQR(A.M12) > t)
                    f = -A.M02 / A.M12;
                norm = 1.0 / Math.Sqrt(1 + SQR(f));
                Q.M00 = norm;
                Q.M10 = f * norm;
                Q.M20 = 0.0;
            }
            else                      // This is the standard branch
            {
                norm = Math.Sqrt(1.0 / norm);
                Q.M00 = Q.M00 * norm;
                Q.M10 = Q.M10 * norm;
                Q.M20 = Q.M20 * norm;
            }


            // Prepare calculation of second eigenvector
            t = w.X - w.Y;
            if (Math.Abs(t) > 8.0 * DBL_EPSILON * wmax)
            {
                // For non-degenerate eigenvalue, calculate second eigenvector by the formula
                //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
                A.M00 += t;
                A.M11 += t;
                Q.M01 = Q.M01 + A.M02 * w.Y;
                Q.M11 = Q.M11 + A.M12 * w.Y;
                Q.M21 = A.M00 * A.M11 - Q.M21;
                norm = SQR(Q.M01) + SQR(Q.M11) + SQR(Q.M21);
                n0 = n0tmp + SQR(A.M00);
                n1 = n1tmp + SQR(A.M11);
                error = n0 * n1;

                if (n0 <= thresh)       // If the first column is zero, then (1,0,0) is an eigenvector
                {
                    Q.M01 = 1.0;
                    Q.M11 = 0.0;
                    Q.M21 = 0.0;
                }
                else if (n1 <= thresh)  // If the second column is zero, then (0,1,0) is an eigenvector
                {
                    Q.M01 = 0.0;
                    Q.M11 = 1.0;
                    Q.M21 = 0.0;
                }
                else if (norm < SQR(64.0 * DBL_EPSILON) * error)
                {                       // If angle between A[0] and A[1] is too small, don't use
                    t = SQR(A.M01);     // cross product, but calculate v ~ (1, -A0/A1, 0)
                    f = -A.M00 / A.M01;
                    if (SQR(A.M11) > t)
                    {
                        t = SQR(A.M11);
                        f = -A.M01 / A.M11;
                    }
                    if (SQR(A.M12) > t)
                        f = -A.M02 / A.M12;
                    norm = 1.0 / Math.Sqrt(1 + SQR(f));
                    Q.M01 = norm;
                    Q.M11 = f * norm;
                    Q.M21 = 0.0;
                }
                else
                {
                    norm = Math.Sqrt(1.0 / norm);
                    Q.M01 = Q.M01 * norm;
                    Q.M11 = Q.M11 * norm;
                    Q.M21 = Q.M21 * norm;
                }
            }
            else
            {
                // For degenerate eigenvalue, calculate second eigenvector according to
                //   v[1] = v[0] x (A - w[1]).e[i]
                //   
                // This would really get to complicated if we could not assume all of A to
                // contain meaningful values.
                A.M10 = A.M01;
                A.M20 = A.M02;
                A.M21 = A.M12;
                A.M00 += w.X;
                A.M11 += w.X;
                for (i = 0; i < 3; i++)
                {
                    A[i, i] -= w.Y;
                    n0 = SQR(A[0, i]) + SQR(A[1, i]) + SQR(A[2, i]);
                    if (n0 > thresh)
                    {
                        Q.M01 = Q.M10 * A[2, i] - Q.M20 * A[1, i];
                        Q.M11 = Q.M20 * A[0, i] - Q.M00 * A[2, i];
                        Q.M21 = Q.M00 * A[1, i] - Q.M10 * A[0, i];
                        norm = SQR(Q.M01) + SQR(Q.M11) + SQR(Q.M21);
                        if (norm > SQR(256.0 * DBL_EPSILON) * n0) // Accept cross product only if the angle between
                        {                                         // the two vectors was not too small
                            norm = Math.Sqrt(1.0 / norm);
                            Q.M01 = Q.M01 * norm;
                            Q.M11 = Q.M11 * norm;
                            Q.M21 = Q.M21 * norm;
                            break;
                        }
                    }
                }

                if (i == 3)    // This means that any vector orthogonal to v[0] is an EV.
                {
                    for (j = 0; j < 3; j++)
                        if (Q[j, 0] != 0.0)                                   // Find nonzero element of v[0] ...
                        {                                                     // ... and swap it with the next one
                            norm = 1.0 / Math.Sqrt(SQR(Q[j, 0]) + SQR(Q[(j + 1) % 3, 0]));
                            Q[j, 1] = Q[(j + 1) % 3, 0] * norm;
                            Q[(j + 1) % 3, 1] = -Q[j, 0] * norm;
                            Q[(j + 2) % 3, 1] = 0.0;
                            break;
                        }
                }
            }

            // Calculate third eigenvector according to
            //   v[2] = v[0] x v[1]
            Q.M02 = Q.M10 * Q.M21 - Q.M20 * Q.M11;
            Q.M12 = Q.M20 * Q.M01 - Q.M00 * Q.M21;
            Q.M22 = Q.M00 * Q.M11 - Q.M10 * Q.M01;
        }

        /// <summary>
        /// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
        /// (unitary) Householder transformations:
        ///           [ d[0]  e[0]       ]
        ///   A = Q . [ e[0]  d[1]  e[1] ] . Q^T
        ///           [       e[1]  d[2] ]
        /// The function accesses only the diagonal and upper triangular
        /// parts of A. The access is read-only.
        /// </summary>
        public static void Dsytrd3(this in M33d A, out M33d Q, out V3d d, out V2d e)
        {
            V3d u, q;
            double omega, f;
            double K, h, g;
            double ex, ey;

            // Initialize Q to the identity matrix.
            Q = M33d.Identity;

            // Bring first row and column to the desired form .
            h = A.M01 * A.M01 + A.M02 * A.M02;
            g = (A.M01 > 0) ? -Math.Sqrt(h) : Math.Sqrt(h);
            ex = g;
            f = g * A.M01;
            u.Y = A.M01 - g;
            u.Z = A.M02;

            omega = h - f;
            if (omega > 0.0)
            {
                omega = 1.0 / omega;
                K = 0.0;
                f = A.M11 * u.Y + A.M12 * u.Z; q.Y = omega * f; K += u.Y * f;
                f = A.M12 * u.Y + A.M22 * u.Z; q.Z = omega * f; K += u.Z * f;

                K *= 0.5 * omega * omega;
                q.Y = q.Y - K * u.Y;
                q.Z = q.Z - K * u.Z;

                d = new V3d(A.M00, A.M11 - 2.0 * q.Y * u.Y, A.M22 - 2.0 * q.Z * u.Z);

                // Store inverse Householder transformation in Q.
                f = omega * u.Y;
                Q.M11 = Q.M11 - f * u.Y;
                Q.M21 = Q.M21 - f * u.Z;
                f = omega * u.Z;
                Q.M12 = Q.M12 - f * u.Y;
                Q.M22 = Q.M22 - f * u.Z;

                // Calculate updated A[1][2] and store it in e[1].
                ey = A.M12 - q.Y * u.Z - u.Y * q.Z;
            }
            else
            {
                d = new V3d(A.M00, A.M11, A.M22);
                ey = A.M12;
            }

            e = new V2d(ex, ey);
        }

        /// <summary>
        /// Finds the three roots w_j of the secular equation
        ///   f(w_j) = 1 + Sum[ z_i / (d_i - w_j) ]  ==  0.
        /// It is assumed that d_0 leq d_1 leq d_2, and that all z_i have the same sign.
        /// The arrays P_i will contain the information required for the calculation
        /// of the eigenvectors:
        ///   P_ij = d_i - w_j.
        /// These differences can be obtained with better accuracy from intermediate results.
        /// </summary>
        public static void Slvsec3(in V3d d, in V3d z, out V3d w, out M33d R, int i0, int i1, int i2)
        {
            V4d a;  // Bounds of the intervals bracketing the roots
            double delta;           // Shift of the d_i which ensures better accuracy
            var dd = V3d.Zero; // Shifted coefficients dd_i = d_i - delta
            double xl, xh;          // Interval which straddles the current root. f(xl) < 0, f(xh) > 0
            double x;               // Current estimates for the root
            V3d x0;      // Analytically calculated roots, used as starting values
            double F, dF;           // Function value f(x) and derivative f'(x)
            double dx, dxold;       // Current and last stepsizes
            double error;           // Numerical error estimate, used for termination condition
            V3d t;                  // Temporary storage used for evaluating f
            double alpha, beta, gamma;       // Coefficients of polynomial f(x) * Product [ d_i - x ]
            double p, sqrt_p, q, c, s, phi;  // Intermediate results of analytical calculation
            w = V3d.Zero;

            // Determine intervals which must contain the roots
            if (z[0] > 0)
            {
                a = new V4d(
                    d[i0], d[i1], d[i2],
                    Math.Abs(d.X + 3.0 * z.X) + Math.Abs(d.Y + 3.0 * z.Y) + Math.Abs(d.Z + 3.0 * z.Z)
                    );
            }
            else
            {
                a = new V4d(
                    -Math.Abs(d.X + 3.0 * z.X) - Math.Abs(d.Y + 3.0 * z.Y) - Math.Abs(d.Z + 3.0 * z.Z),
                    d[i0], d[i1], d[i2]
                    );
            }

            // Calculate roots of f(x) = 0 analytically (analogous to ZHEEVC3)
            t = new V3d(d.Y * d.Z, d.X * d.Z, d.X * d.Y);
            gamma = t.X * d.X + (z.X * t.X + z.Y * t.Y + z.Z * t.Z);    // Coefficients
            beta = (z.X * (d.Y + d.Z) + z.Y * (d.X + d.Z) + z.Z * (d.X + d.Y))
                     + (t.X + t.Y + t.Z);
            alpha = (z.X + z.Y + z.Z) + (d.X + d.Y + d.Z);

            p = alpha * alpha - 3.0 * beta;    // Transformation that removes the x^2 term
            q = alpha * (p - (3.0 / 2.0) * beta) + (27.0 / 2.0) * gamma;
            sqrt_p = Math.Sqrt(Math.Abs(p));

            phi = 27.0 * (0.25 * beta * beta * (p - beta) - gamma * (q - 27.0 / 4.0 * gamma));
            phi = (1.0 / 3.0) * Math.Atan2(Math.Sqrt(Math.Abs(phi)), q);
            c = sqrt_p * Math.Cos(phi);
            s = (1.0 / M_SQRT3) * sqrt_p * Math.Abs(Math.Sin(phi));

            x0 = new V3d((1.0 / 3.0) * (alpha - c));
            if (c > s)  // Make sure the roots are in ascending order.
            {
                x0.X -= s;
                x0.Y += s;
                x0.Z += c;
            }
            else if (c < -s)
            {
                x0.X += c;
                x0.Y -= s;
                x0.Z += s;
            }
            else
            {
                x0.X -= s;
                x0.Y += c;
                x0.Z += s;
            }

            // Refine roots with a combined Bisection/Newton-Raphson method.
            R = M33d.Zero;
            for (int i = 0; i < 3; i++)
            {
                xl = a[i];          // Lower bound of bracketing interval.
                xh = a[i + 1];      // Upper bound of bracketing interval.
                dx = dxold = 0.5 * (xh - xl);

                // Make sure that xl != xh
                if (dx == 0.0)
                {
                    w[i] = xl;
                    for (int j = 0; j < 3; j++)
                        R[j, i] = d[j] - xl;
                    continue;
                }

                // Shift the root close to zero to achieve better accuracy.
                if (x0[i] >= xh)
                {
                    delta = xh;
                    x = -dx;
                    for (int j = 0; j < 3; j++)
                    {
                        dd[j] = d[j] - delta;
                        R[j, i] = dd[j] - x;
                    }
                }
                else if (x0[i] <= xl)
                {
                    delta = xl;
                    x = dx;
                    for (int j = 0; j < 3; j++)
                    {
                        dd[j] = d[j] - delta;
                        R[j, i] = dd[j] - x;
                    }
                }
                else
                {
                    delta = x0[i];
                    x = 0.0;
                    for (int j = 0; j < 3; j++)
                        R[j, i] = dd[j] = d[j] - delta;
                }
                xl -= delta;
                xh -= delta;

                // Make sure that f(xl) < 0 and f(xh) > 0 .
                if (z.X < 0.0)
                {
                    var tmp = xh;
                    xh = xl;
                    xl = tmp;
                }

                // Main iteration loop
                for (int nIter = 0; nIter < 500; nIter++)
                {
                    // Evaluate f and f', and calculate an error estimate.
                    F = 1.0;
                    dF = 0.0;
                    error = 1.0;
                    for (int j = 0; j < 3; j++)
                    {
                        t.X = 1.0 / R[j, i];
                        t.Y = z[j] * t[0];
                        t.Z = t[1] * t[0];
                        F += t[1];
                        error += Math.Abs(t[1]);
                        dF += t[2];
                    }

                    // Check for convergence 
                    if (Math.Abs(F) <= DBL_EPSILON * (8.0 * error + Math.Abs(x * dF)))
                        break;

                    // Adjust interval boundaries
                    if (F < 0.0)
                        xl = x;
                    else
                        xh = x;

                    // Check, whether Newton-Raphson would converge fast enough. If so,
                    // give it a try. If not, or if it would run out of bounds, use bisection
                    if (Math.Abs(2.0 * F) < Math.Abs(dxold * dF))
                    {
                        dxold = dx;
                        dx = F / dF;
                        x = x - dx;
                        if ((x - xh) * (x - xl) >= 0.0)
                        {
                            dx = 0.5 * (xh - xl);
                            x = xl + dx;
                        }
                    }
                    else
                    {
                        dx = 0.5 * (xh - xl);
                        x = xl + dx;
                    }

                    // Prepare next iteration
                    for (int j = 0; j < 3; j++)
                        R[j, i] = dd[j] - x;
                }

                // Un-shift result
                w[i] = x + delta;
            }
        }







        
        private static void Dsyev2(
            in double A, in double B, in double C,
            out double rt1, out double rt2, out double cs, out double sn
            )
        {
            double sm = A + C;
            double df = A - C;
            double rt = Math.Sqrt(df * df + 4.0 * B * B);
            double t;

            if (sm > 0.0)
            {
                rt1 = 0.5 * (sm + rt);
                t = 1.0 / rt1;
                rt2 = (A * t) * C - (B * t) * B;
            }
            else if (sm < 0.0)
            {
                rt2 = 0.5 * (sm - rt);
                t = 1.0 / rt2;
                rt1 = (A * t) * C - (B * t) * B;
            }
            else // This case needs to be treated separately to avoid div by 0.
            {
                rt1 = 0.5 * rt;
                rt2 = -0.5 * rt;
            }

            // Calculate eigenvectors.
            cs = (df > 0.0) ? (df + rt) : (df - rt);

            if (Math.Abs(cs) > 2.0 * Math.Abs(B))
            {
                t = -2.0 * B / cs;
                sn = 1.0 / Math.Sqrt(1.0 + t * t);
                cs = t * sn;
            }
            else if (Math.Abs(B) == 0.0)
            {
                cs = 1.0;
                sn = 0.0;
            }
            else
            {
                t = -0.5 * cs / B;
                cs = 1.0 / Math.Sqrt(1.0 + t * t);
                sn = t * cs;
            }

            if (df > 0.0)
            {
                t = cs;
                cs = -sn;
                sn = t;
            }
        }
        
        private const double DBL_EPSILON = 2.2204460492503131e-16;
        private const float FLT_EPSILON = 1.19209290e-07F;
        private const double M_SQRT3 = 1.73205080756887729352744634151;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double SQR(double x) => x * x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double[,] Unpack(M33d m) => new[,] { { m.M00, m.M01, m.M02 }, { m.M10, m.M11, m.M12 }, { m.M20, m.M21, m.M22 } };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double[] Unpack(V3d v) => new[] { v.X, v.Y, v.Z };
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double[] Unpack(V2d v) => new[] { v.X, v.Y };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static M33d Pack(double[,] m) => new M33d(m[0, 0], m[0, 1], m[0, 2], m[1, 0], m[1, 1], m[1, 2], m[2, 0], m[2, 1], m[2, 2]);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static V3d Pack(double[] v) => new V3d(v[0], v[1], v[2]);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static V2d PackV2d(double[] v) => new V2d(v[0], v[1]);
    }
}
