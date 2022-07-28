namespace CubicSplineInterpolation;

public class Point
{
    public Point(double x, double y)
    {
        X = x;
        Y = y;
    }

    public double X { get; }
    public double Y { get; }
}

public class F
{
    public F(decimal a, decimal b, decimal c, decimal d, FRange range)
    {
        A = a;
        B = b;
        C = c;
        D = d;
        Range = range;
    }

    public class FRange
    {
        public FRange(decimal xMin, decimal xMax)
        {
            XMin = xMin;
            XMax = xMax;
        }

        public decimal XMin { get; }
        public decimal XMax { get; }
    }

    public decimal A { get; }
    public decimal B { get; }
    public decimal C { get; }
    public decimal D { get; }
    public FRange Range { get; }
}

public enum Boundary
{
    Quadratic,
    Notaknot,
    Periodic,
    Natural
}

public static class Helpers
{
    public static IEnumerable<F> Interpolate(IEnumerable<Point> points, Boundary boundary)
    {
        int row = 0;
        int solutionIndex = (points.Count() - 1) * 4;

        // initialize matrix
        var m = new decimal[solutionIndex][]; // rows
        for (var i = 0; i < solutionIndex; i++)
        {
            // columns (rows + 1)
            m[i] = new decimal[solutionIndex + 1];
        }

        // splines through p equations
        for (var functionNr = 0; functionNr < points.Count() - 1; functionNr++, row++)
        {
            var p0 = points.ElementAt(functionNr);
            var p1 = points.ElementAt(functionNr + 1);
            m[row][functionNr * 4 + 0] = (decimal)Math.Pow(p0.X, 3);
            m[row][functionNr * 4 + 1] = (decimal)Math.Pow(p0.X, 2);
            m[row][functionNr * 4 + 2] = (decimal)p0.X;
            m[row][functionNr * 4 + 3] = 1;
            m[row][solutionIndex] = (decimal)p0.Y;

            m[++row][functionNr * 4 + 0] = (decimal)Math.Pow(p1.X, 3);
            m[row][functionNr * 4 + 1] = (decimal)Math.Pow(p1.X, 2);
            m[row][functionNr * 4 + 2] = (decimal)p1.X;
            m[row][functionNr * 4 + 3] = 1;
            m[row][solutionIndex] = (decimal)p1.Y;
        }

        // first derivative
        for (var functionNr = 0; functionNr < points.Count() - 2; functionNr++, row++)
        {
            var p1 = points.ElementAt(functionNr + 1);
            m[row][functionNr * 4 + 0] = 3 * (decimal)Math.Pow(p1.X, 2);
            m[row][functionNr * 4 + 1] = 2 * (decimal)p1.X;
            m[row][functionNr * 4 + 2] = 1;
            m[row][functionNr * 4 + 4] = -3 *(decimal)Math.Pow(p1.X, 2);
            m[row][functionNr * 4 + 5] = -2 * (decimal)p1.X;
            m[row][functionNr * 4 + 6] = -1;
        }

        // second derivative
        for (var functionNr = 0; functionNr < points.Count() - 2; functionNr++, row++)
        {
            var p1 = points.ElementAt(functionNr + 1);
            m[row][functionNr * 4 + 0] = 6 * (decimal)p1.X;
            m[row][functionNr * 4 + 1] = 2;
            m[row][functionNr * 4 + 4] = -6 * (decimal)p1.X;
            m[row][functionNr * 4 + 5] = -2;
        }

        // boundary conditions
        switch (boundary)
        {
            case Boundary.Quadratic: // first and last spline quadratic
                m[row++][0] = 1;
                m[row++][solutionIndex - 4 + 0] = 1;
                break;

            case Boundary.Notaknot: // Not-a-knot spline
                m[row][0 + 0] = 1;
                m[row++][0 + 4] = -1;
                m[row][solutionIndex - 8 + 0] = 1;
                m[row][solutionIndex - 4 + 0] = -1;
                break;

            case Boundary.Periodic: // periodic function
                // first derivative of first and last point equal
                m[row][0] = 3 * (decimal)Math.Pow(points.ElementAt(0).X, 2);
                m[row][1] = 2 * (decimal)points.ElementAt(0).X;
                m[row][2] = 1;
                m[row][solutionIndex - 4 + 0] = -3 * (decimal)Math.Pow(points.ElementAt(points.Count() - 1).X, 2);
                m[row][solutionIndex - 4 + 1] = -2 * (decimal)points.ElementAt(points.Count() - 1).X;
                m[row++][solutionIndex - 4 + 2] = -1;

                // second derivative of first and last point equal
                m[row][0] = 6 * (decimal)points.ElementAt(0).X;
                m[row][1] = 2;
                m[row][solutionIndex - 4 + 0] = -6 * (decimal)points.ElementAt(points.Count() - 1).X;
                m[row][solutionIndex - 4 + 1] = -2;
                break;

            case Boundary.Natural: // natural spline
                m[row][0 + 0] = 6 * (decimal)points.ElementAt(0).X;
                m[row++][0 + 1] = 2;
                m[row][solutionIndex - 4 + 0] = 6 * (decimal)points.ElementAt(points.Count() - 1).X;
                m[row][solutionIndex - 4 + 1] = 2;
                break;

            default: throw new Exception("Unknown boundary " + boundary);
        }


        var reducedRowEchelonForm = ReduceRowEchelonForm(m);
        var coefficients = new List<decimal>();

        for (var i = 0; i < reducedRowEchelonForm.Count(); i++)
        {
            coefficients.Add(reducedRowEchelonForm[i][reducedRowEchelonForm[i].Count() - 1]);
        }

        var functions = new List<F>();

        for (var i = 0; i < coefficients.Count(); i += 4)
        {
            functions.Add(
                new F
                (
                    a: coefficients[i],
                    b: coefficients[i + 1],
                    c: coefficients[i + 2],
                    d: coefficients[i + 3],
                    range: new F.FRange((decimal)points.ElementAt(i / 4).X, (decimal)points.ElementAt(i / 4 + 1).X)
                )
            );
        }

        return functions;
    }

    // https://rosettacode.org/wiki/Reduced_row_echelon_form
    public static decimal[][] ReduceRowEchelonForm(decimal[][] mat)
    {
        var lead = 0;
        for (var r = 0; r < mat.Count(); r++)
        {
            if (mat[0].Count() <= lead)
            {
                return mat;
            }
            var i = r;
            while (mat[i][lead] == 0)
            {
                i++;
                if (mat.Count() == i)
                {
                    i = r;
                    lead++;
                    if (mat[0].Count() == lead)
                    {
                        return mat;
                    }
                }
            }

            var tmp = mat[i];
            mat[i] = mat[r];
            mat[r] = tmp;

            var val = mat[r][lead];
            for (var j = 0; j < mat[0].Count(); j++)
            {
                mat[r][j] = mat[r][j] / val;
            }

            for (i = 0; i < mat.Count(); i++)
            {
                if (i == r) continue;
                val = mat[i][lead];
                for (var j = 0; j < mat[0].Count(); j++)
                {
                    mat[i][j] = mat[i][j] - val * mat[r][j];
                }
            }
            lead++;
        }
        return mat;
    }

    public static decimal? CalculateY(IEnumerable<F> fns, decimal x)
    {
        for (var i = 0; i < fns.Count(); i++)
        {
            var f = fns.ElementAt(i);
            if (f.Range.XMin <= x && f.Range.XMax >= x)
            {
                return f.A * x * x * x + f.B * x * x + f.C * x + f.D;
            }
        }
        return null;
    }
}