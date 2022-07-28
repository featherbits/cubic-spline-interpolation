using CubicSplineInterpolation;

var points = new Point[]
{
    new Point(1.2695, 10),
    new Point(1.4060, 15),
    new Point(1.7100, 30),
    new Point(2.1563, 60),
    new Point(2.7522, 120),
    new Point(3.5070, 250),
    new Point(3.9393, 370),
    new Point(4.2349, 480),
    new Point(4.5669, 640),
};

var fns = Helpers.Interpolate(points, Boundary.Natural);

Console.WriteLine(Helpers.CalculateY(fns, 4.6m));