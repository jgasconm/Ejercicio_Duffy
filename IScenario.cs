using System;
using System.Collections.Generic;
using System.Text;
using System.Data.OleDb;
using System.Data;

/// prueba de commit con tortoisegit
namespace DuffyExercise
{
    using ExtensionsRandom;

    public class ScenarioKey //Key class to look for a value into the scenario Dictionary
    {
        public string _ID;
        public double _date;
        public int _path;

        public bool Equals(ScenarioKey Key)
        {
            return (Key._ID == _ID && Key._date == _date && Key._path == _path);
        }

        public override bool Equals(Object obj)
        {
            return Equals(obj as ScenarioKey);
        }
        public override int GetHashCode()
        {
            return _ID.GetHashCode() ^ _date.GetHashCode() ^ _path.GetHashCode();
        }

        public ScenarioKey(string ID, double date, int path)
        {
            _ID = ID;
            _date = date;
            _path = path;
        }

    }

    // BROWNIAN MANAGER ---------------------------------------------------------------------------
    /// <summary>
    ///     Brownian Dispatcher. Given certain logic it delivers brownian motion to every risk Factor.
    ///     Implement the dispatch method.
    /// </summary>
    public class BrownianManager
    {
        Dictionary<string, int> NumFactorsForEachID;
        IGaussianGenerator _gaussianGenerator;

        public BrownianManager(IGaussianGenerator gaussianGenerator)
        {
            _gaussianGenerator = gaussianGenerator;
            NumFactorsForEachID = new Dictionary<string, int>();
        }

        /// <summary>
        ///     Implements the logic for dispatching.
        /// </summary>
        public void initialize(List<RiskFactor> riskFactors)
        {
            foreach(RiskFactor rf in riskFactors)
            {
                NumFactorsForEachID.Add(rf.ID, rf.size());
            }
        }

        public double[] dispatch(string sID)
        {
            return _gaussianGenerator.nextSequence(NumFactorsForEachID[sID]);
        }
    }

    // --------------------------------------------------------------------------------------


    // SCENARIO ---------------------------------------------------------------------------
    /// <summary>
    /// Interface for the Scenario Container. Must provide accesors to the economy data given a filtration.
    /// </summary>
    public interface IScenario
    {
        /// <summary>
        ///     Returns a value for a scalar risk factor for date and path.
        /// </summary>
        double getValue(string ID, double date, int path);
        void setValue(string ID, double date, int path, double value);
    }

    /// <summary>
    ///     Class implementing the Scenario Interface
    /// </summary>
    public class Scenario : IScenario
    {
        //->ATTRIBUTES
        //private List< Dictionary<string, Dictionary<double, double>> >_scenarios;
        private Dictionary<ScenarioKey, double> _scenarios;

        #region IScenario Members

        public double getValue(string ID, double date, int path)
        {
            try
            {
                return _scenarios[new ScenarioKey(ID, date, path)];
            }
            catch
            {
                throw new Exception("Scenario not found.");
            }

        }

        public Scenario()
        {
            /*List<double> Dates = new List<double>(new double[] { 40000, 40100, 40200, 40300, 40400, 40500, 40600, 40700, 40800, 40900 });
            double DFstep = 0.015;            
            _scenarios = new Dictionary<ScenarioKey, double>();
            //IR Curve

            for (int i = 0; i < Dates.Count; i++)
                _scenarios.Add(new ScenarioKey("IRCurve", Dates[i], 0), 1 - i * DFstep);*/
        }

        public void setValue(string ID, double date, int path, double value)
        {

            _scenarios.Add(new ScenarioKey(ID, date, path), value);
        }

        #endregion
    }

    // --------------------------------------------------------------------------------------



    // RISK FACTORS ---------------------------------------------------------------------------
    /// <summary>
    /// Wrapper to a SDE. Encapsulates the dynamic of a Risk Factor 
    /// </summary>
    public abstract class RiskFactor
    {
        // -> ATTRIBUTES
        // --------------------
        private string _ID;

        // -> CONSTRCTOR
        // --------------------
        protected RiskFactor(string ID)
        {
            _ID = ID;
        }
        protected RiskFactor() {}

        // -> ACCESORS
        // ------------------
        public string ID { get { return _ID; } }


        // -> ABSTRACT METHODS
        //---------------------------------------
        public void evolve(double dateFrom, double dateTo, int path, IScenario scenario,
                                    double[] correlatedBrownians, BrownianManager brownianDispacther)
        {
            evolve(dateFrom, dateTo, path, scenario, brownianDispacther.dispatch(ID));
        }
        public abstract void evolve(double dateFrom, double dateTo, int path,
                                    IScenario scenario, double[] correlatedBrownians);
        /// <summary>
        ///     Returns the number of Brownian Motions needed to evolve the Risk Factor
        /// </summary>
        public abstract int size();
    }


    /// <summary>
    ///     Equity Dynamics: Implement Contructor, Abstract Method ...
    /// </summary>
    public class EquityRiskFactor : RiskFactor
    {
        //ATTRIBUTES
        private double _vol;
        private double _riskFreeRate;

        //CONSTRUCTOR
        public EquityRiskFactor(string ID) : base(ID)
        {
            throw new Exception("EquityRiskFactor constructor must be called including values for vol and riskFreeRate!");
        }

        public EquityRiskFactor(string ID, double vol, double riskFreeRate) : base(ID)
        {
            _vol= vol;
            _riskFreeRate = riskFreeRate;
        }

        // -> METHODS
        // --------------------
        public override void evolve(double dateFrom, double dateTo, int path,
                                    IScenario scenario, double[] correlatedBrownians)
        {
            double currentValue;
            double nextValue;
            double incT = (dateTo - dateFrom) / 365.25;
            try
            {                
                currentValue = scenario.getValue(ID, dateFrom, path);
                nextValue = currentValue * Math.Exp(_riskFreeRate*incT) 
                                         * Math.Exp(- (Math.Pow(_vol, 2) * (incT/2)) )
                                         * Math.Exp(_vol * Math.Sqrt(incT) * correlatedBrownians[0]);
                scenario.setValue(ID, dateTo, path, nextValue);
            }
            catch
            {
                throw new Exception("Error in method evolve of risk factor " + ID + " from date " + dateFrom);
            }
        }

        public override int size()
        {
            return 1;
        }
    }

    // --------------------------------------------------------------------------------------


    // GAUSSIAN GENERATOR ---------------------------------------------------------------------------

    /// <summary>
    ///  Gaussian generator ...
    /// </summary>
    public interface IGaussianGenerator
    {
        double next();
        double[] nextSequence(int dim);
    }

    namespace ExtensionsRandom
    {
        /// <summary>
        ///  Extending class Random (with "cow" method NextDouble) to have method NextGaussian
        /// </summary>
        public static class RandomExtensions
        {
            public static double NextGaussian(this System.Random rand, double mu = 0.0, double sigma = 1.0)
            {
                double u1 = rand.NextDouble();
                double u2 = rand.NextDouble();

                double rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
                double rand_normal = mu + sigma * rand_std_normal;

                return rand_normal;
            }
        }
    }

    public class BoxMuller : IGaussianGenerator
    {
        public double next()
        {
            Random rand = new System.Random();
            return rand.NextGaussian();
        }

        public double[] nextSequence(int dim)
        {
            double[] sequence = new double[dim];

            for (int j = 0; j < dim; ++j)
            {
                sequence.SetValue(next(), j);
            }

            return sequence;
        }
    }

    // --------------------------------------------------------------------------------------


    // MONTE CARLO ENGINE -------------------------------------------------------------------
    public abstract class Correlator
    {
        // This method has as input independent normals, correlates them and then multiplies by sqrt(deltaT) to have the Brownian Motion in the interval
        public abstract double[] correlate(double dateFrom, double dateTo, double[] independentBM);
    }

    public class Cholesky : Correlator
    {
        double[,] _correlMatrix;

        public double[,] CorrelMatrix
        {
            get {return _correlMatrix;}
            set {_correlMatrix = value;}
        }
        public Cholesky(double[,] correlMatrix)
        {
            _correlMatrix = correlMatrix;
        }


        public override double[] correlate(double dateFrom, double dateTo, double[] independentBM)
        {
            double[] correlatedBrownians;

            ///
            /// *^*�* MISSING CHOLESKY DECOMPOSITION AND BASE CHANGE PROCESS
            /// 

            double sqrtDeltaT = Math.Sqrt( (dateTo-dateFrom) / 365.25 );

            correlatedBrownians = choleskyMultiplication(independentBM);
            for (int i = 0; i < correlatedBrownians.Length; i++)
            {
                correlatedBrownians[i] *= sqrtDeltaT;
            }

            return correlatedBrownians;
        }

        public double[] choleskyMultiplication(double[] independentNormals)
        {
            ///
            /// *^*�* MISSING CHOLESKY DECOMPOSITION AND BASE CHANGE PROCESS
            /// 
            return independentNormals; 
        }
    }

    /// <summary>
    ///     Class that evolves every Risk Factor .... 
    /// </summary>
    public class Evolver
    {
        // -> ATTRIBUTES
        // ---------------------
        List<RiskFactor> _l_RiskFactors;
        IGaussianGenerator _gaussianGenerator;
        BrownianManager _brownianManager;
        Correlator _correlator;

        //CONSTRUCTOR

        public Evolver(List<RiskFactor> Equities, double [,] correlMatrix)
        {
            _gaussianGenerator = new BoxMuller();
            _brownianManager = new BrownianManager(_gaussianGenerator);
            _correlator = new Cholesky(correlMatrix);
            _l_RiskFactors = Equities;
            _brownianManager.initialize(_l_RiskFactors);
        }
        // -> METHODS
        // ---------------------

        /// <summary>
        ///     Implement: evolve every Risk Factor from dateFrom to dateTo.
        /// </summary>
        public void evolve(double dateFrom, double dateTo, int path, IScenario scenario)
        {
            double[] independentNormals; // Just the independent gaussian random numer
            double[] correlatedBrownians; // Correlated brownian motion, including sqrt(deltaT)

            foreach(RiskFactor Factor in _l_RiskFactors)
            {
                independentNormals = _brownianManager.dispatch(Factor.ID);
                correlatedBrownians = _correlator.correlate(dateFrom, dateTo, independentNormals);
                Factor.evolve(dateFrom, dateTo, path, scenario, correlatedBrownians);
            }
        }
    }

    // --------------------------------------------------------------------------------------



    // INSTRUMENT AND PRICER ---------------------------------------------------------------------------

    /// <summary>
    ///     Encapsulates both the contractual info and the pricer.
    ///     It offers the price interface.
    /// </summary>
    public class Instrument
    {
        // -> ATTRIBUTES
        // -------------------
        Pricer _pricer;
        InstrumentContract _InstrumentContract;

        // -> METHODS 
        // --------------------
        public virtual double price(double date, int path, IScenario scenario)
        {
            return _pricer.price(date, path, _InstrumentContract, scenario);
        }


        // -> ACCEORS & MODIFIERS
        // -----------------------
        public Pricer pricer { get { return _pricer; } set { _pricer = value; } }
        public InstrumentContract instrumentContract { get { return _InstrumentContract; } set { _InstrumentContract = value; } }
    }
    /// <summary>
    ///     Information about the contractual agreement.
    ///     Populate with NEEDED ATTRIBUTES
    /// </summary>
    public abstract class InstrumentContract
    {
        // -> ATTRIBUTES
        // -------------------

    }
    /// <summary>
    ///     Implementst the logic of pricing.
    /// </summary>
    public abstract class Pricer
    {
        public abstract double price(double date, int path, InstrumentContract instrumentContract, IScenario scenario);
    }

    // --------------------------------------------------------------------------------------




    // METRICS (HISTOGRAM, ETC ..)  -----------------------------------------------------------------------------------

    public abstract class Metric
    {
        public abstract void addNPVToMetric(double date, int path, double NPV, IScenario scenario);
        public abstract void writeToConsole();
    }
    /// <summary>
    ///     Stores the Distribution for the NPV for a certain date
    /// </summary>
    public class Histogram : Metric
    {
        public Histogram(List<double> distribution)
        {

        }

        public override void addNPVToMetric(double date, int path, double NPV, IScenario scenario)
        {
            throw new Exception("The method or operation is not implemented.");
        }
        public override void writeToConsole()
        {
            throw new Exception("The method or operation is not implemented.");
        }
    }

    // --------------------------------------------------------------------------------------



    // AGGREGATOR -----------------------------------------------------------------------------------

    /// <summary>
    ///     Main Object: Creates the scenario, prices every instrument according to that scenario and stores the NPV's distribution
    ///         Must be capable to return percentiles for every date.
    /// </summary>
    public class MCEngine
    {
        // -> ATTRIBUTES
        // --------------------
        Evolver _engine;
        List<Instrument> _l_Instruments;
        List<Metric> _l_Metric;

        // CONSTRUCTOR (TO BE IMPLEMENTED)
        // -------------------------------
        public MCEngine()
        {
            Initialize();
            _l_Instruments = new List<Instrument>();
            _l_Metric = new List<Metric>();

        }

        public void Initialize() //Read Excel info 
        {
            List<double> SimDates;
            List<RiskFactor> Equities;
            double[,] CorrelationMatrix;
            SimDates = ReadSimulationDates();
            Equities = ReadEquities();
            CorrelationMatrix = ReadCorrelationMatrix();
            _engine = new Evolver(Equities, CorrelationMatrix);


        }

        public List<double> ReadSimulationDates()
        {
            List<double> Dates = new List<double>();
            ReadRange RangeHandle = new ReadRange();
            DataSet Dates_Data = RangeHandle.Read("SimDates");
            DataRowCollection Row = Dates_Data.Tables["SimDates"].Rows;

            foreach(DataRow DRow in Row)
            {
                Dates.Add(Convert.ToDouble(DRow[0]));
            }

            return Dates;

        }

        public List<RiskFactor> ReadEquities()
        {
            List<RiskFactor> Equities = new List<RiskFactor>();
            ReadRange RangeHandle = new ReadRange();
            DataSet Equity_Data = RangeHandle.Read("EquityFeatures");
            DataSet Rate_Data = RangeHandle.Read("Rates");
            DataRowCollection EquityRows = Equity_Data.Tables["EquityFeatures"].Rows;
            DataColumnCollection EquityColumns = Equity_Data.Tables["EquityFeatures"].Columns;

            DataRowCollection RateRows = Rate_Data.Tables["Rates"].Rows;

            foreach (DataRow DRow in EquityRows)
            {
                Equities.Add(new EquityRiskFactor(Convert.ToString(DRow[0]), Convert.ToDouble(DRow[1]), Convert.ToDouble(RateRows[0][1])));
            }
            return Equities;

        }

        public double[,] ReadCorrelationMatrix()
        {
            double[,] correlMatrix;
            ReadRange Matrix_Range = new ReadRange();
            DataSet Matrix_Data = Matrix_Range.Read("MatrixCorrelation");
            DataColumnCollection Column = Matrix_Data.Tables["MatrixCorrelation"].Columns;
            DataRowCollection Row = Matrix_Data.Tables["MatrixCorrelation"].Rows;

            correlMatrix = new double[Row.Count, Column.Count];

            for (int i = 0; i < Row.Count; i++)
            {
                for (int j = 0; j < Column.Count; j++)
                {
                    correlMatrix[i, j] = Convert.ToDouble(Row[i][j]);
                }
            }

            return correlMatrix;
        }

        // -> METHODS
        // -------------------------------
        public void calculate(List<double> dates, int noSim)
        {
            Scenario scenario = new Scenario();           

            // -> CREATE SCENARIO ...
            for (int j = 0; j < noSim; ++j)
            {
                double EquitySpot = 9.89; //BBVA
                scenario.setValue("BBVA", dates[0], j, EquitySpot);

                for (int i = 1; i < dates.Count; ++i)
                {
                    // -> EVOLVE ..
                    _engine.evolve(dates[i - 1], dates[i], j, scenario);

                    // -> PRICE
                    double npv = 0.0;
                    for (int k = 0; k < _l_Instruments.Count; ++k)
                        npv += _l_Instruments[k].price(dates[i], j, scenario);

                    // -> CALL TO METRIC TO TAKE INTO ACCOUNT THE NPV
                    foreach (Metric metric in _l_Metric)
                        metric.addNPVToMetric(dates[i], j, npv, scenario);
                }
            }

        }

       
    }

    public class ReadRange
    {
        public DataSet Read(string RangeName)
        {
            String sConnectionString = "Provider=Microsoft.Jet.OLEDB.4.0;Data Source=C:/Users/e024097/Documents/Visual Studio 2013/Projects/ConsoleApplication1/Book1.xls;Extended Properties=Excel 8.0;";
            OleDbConnection objConn = new OleDbConnection(sConnectionString);
            objConn.Open();
            OleDbCommand objCmd = new OleDbCommand("SELECT * FROM " + RangeName, objConn);
            OleDbDataAdapter Adapter = new OleDbDataAdapter(objCmd);
            DataSet Data = new DataSet();
            Adapter.Fill(Data, RangeName);
            objConn.Close();
            return Data;

        }
    }

    // --------------------------------------------------------------------------------------


}
