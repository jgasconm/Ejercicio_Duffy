using System;
using System.Collections.Generic;
using System.Text;


namespace DuffyExercise
{
    public class Exe
    {
        

        public static void Main()
        {
            MCEngine _MCEngine = new MCEngine();
            List<double> Dates = new List<double>(new double[]{40000, 40100, 40200, 40300, 40400, 40500, 40600, 40700, 40800, 40900 });
            int noSim = 10;
            _MCEngine.calculate(Dates, noSim);
        }
    }
}
