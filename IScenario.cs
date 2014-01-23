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

    public class NPVKey //Key class to look for a value into the NPV Dictionary
    {
        public double _date;
        public int _path;

        public bool Equals(NPVKey Key)
        {
            return (Key._date == _date && Key._path == _path);
        }

        public override bool Equals(Object obj)
        {
            return Equals(obj as NPVKey);
        }
        public override int GetHashCode()
        {
            return _date.GetHashCode() ^ _path.GetHashCode();
        }

        public NPVKey(double date, int path)
        {
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

        public double[] dispatch(int numFactors)
        {
            return _gaussianGenerator.nextSequence(numFactors);
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

        public Scenario(List<RiskFactor> RFList, double StartDate)
        {
            _scenarios = new Dictionary<ScenarioKey, double>();
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
        private double _spot;

        // -> CONSTRCTOR
        // --------------------
        protected RiskFactor(string ID, double spot)
        {
            _ID = ID;
            _spot = spot;
        }
        protected RiskFactor() {}

        // -> ACCESORS
        // ------------------
        public string ID { get { return _ID; } }
        public double spot { get { return _spot; } }


        // -> ABSTRACT METHODS
        //---------------------------------------
        /*public void evolve(double dateFrom, double dateTo, int path, IScenario scenario,
                                    double[] correlatedBrownians, BrownianManager brownianDispacther)
        {
            evolve(dateFrom, dateTo, path, scenario, brownianDispacther.dispatch(ID));
        }*/
        public abstract void evolve(double dateFrom, double dateTo, int path,
                                    IScenario scenario, double correlatedBrownians);

        public override bool Equals(object obj)
        {
            return base.Equals(obj as string);
        }

        public bool Equals(string Factor)
        {
            return Factor == ID;
        }    
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

        public double vol{get{return _vol;}}
        public double riskFreeRate{get{return _riskFreeRate;}}

        //CONSTRUCTOR
        public EquityRiskFactor(string ID, double spot) : base(ID, spot)
        {
            throw new Exception("EquityRiskFactor constructor must be called including values for vol and riskFreeRate!");
        }

        public EquityRiskFactor(string ID, double spot, double vol, double riskFreeRate) : base(ID, spot)
        {
            _vol= vol;
            _riskFreeRate = riskFreeRate;
        }

        // -> METHODS
        // --------------------
        public override void evolve(double dateFrom, double dateTo, int path,
                                    IScenario scenario, double correlatedBrownians)
        {
            double currentValue;
            double nextValue;
            double incT = (dateTo - dateFrom) / 365.25;
            try
            {                
                currentValue = scenario.getValue(ID, dateFrom, path);
                nextValue = currentValue * Math.Exp(_riskFreeRate*incT) 
                                         * Math.Exp(- (Math.Pow(_vol, 2) * (incT/2)) )
                                         * Math.Exp(_vol * Math.Sqrt(incT) * correlatedBrownians);
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
        private static Random _rand;

        public BoxMuller()
        {
            _rand = new System.Random();
        }
        public double next()
        {
            return _rand.NextGaussian();
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

        // Input: independent NORMALS
        // Output: correlated BROWNIANS (normal * sqrt(deltaT))
        public override double[] correlate(double dateFrom, double dateTo, double[] independentNormals)
        {
            double[] correlatedBrownians;
            double[][] _baseChangeMatrix;

            _baseChangeMatrix = choleskyFactorization(_correlMatrix);
            correlatedBrownians = baseChangeMultiplication(_baseChangeMatrix, independentNormals);

            double sqrtDeltaT = Math.Sqrt((dateTo - dateFrom) / 365.25);
            
            for (int i = 0; i < correlatedBrownians.Length; i++)
            {
                correlatedBrownians[i] *= sqrtDeltaT;
            }

            return correlatedBrownians;
        }

        // Supposing a square correlation matrix!
        public double[] baseChangeMultiplication(double[][] baseChangeMatrix, double[] independentNormals)
        {
            int dimens = CorrelMatrix.GetLength(0);
            double[] correlatedNormals = new double[dimens];

            for (int i = 0; i < dimens; i++)
            {
                double sum = 0.0;
                for (int j = i; j < dimens; j++)
                {
                    if(j <= i)
                        sum += baseChangeMatrix[i][j] * independentNormals[j];
                    else
                        sum += baseChangeMatrix[j][i] * independentNormals[j];
                }
                correlatedNormals[i] = sum;
            }

            return correlatedNormals;
        }

        public double[][] choleskyFactorization(double[,] CorrelMatrix)
        {
            // Just uses lower triangular correlation matrix.
            // Therefore a triangular jagged array can be received 
            // (although the routine works with a square array).

            // It returns a lower triangular jagged array.

            double sum;
            int k;
            int dimens = CorrelMatrix.GetLength(0);

            double[][] result = new double[dimens][];
            for (int i = 0; i < dimens; i++) result[i] = new double[i + 1];

            for (int i = 0; i < dimens; i++)
            {
                for (int j = i; j < dimens; j++)
                {
                    for (sum = CorrelMatrix[j,i], k = i - 1; k >= 0; k--) sum -= result[i][k] * result[j][k];
                    if (i == j)
                    {
                        if (sum <= 0.0)
                            throw new Exception("Non semi-positive definite matrix.");
                        result[i][i] = Math.Sqrt(sum);
                    }
                    else result[j][i] = sum / result[i][i];
                }
            }

            return result;
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
        public List<RiskFactor> RiskFactors { get { return _l_RiskFactors; } }
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

            
            independentNormals = _brownianManager.dispatch(_l_RiskFactors.Count);
            correlatedBrownians = _correlator.correlate(dateFrom, dateTo, independentNormals);
            for (int i = 0; i < _l_RiskFactors.Count; i++ )
                _l_RiskFactors[i].evolve(dateFrom, dateTo, path, scenario, correlatedBrownians[i]);
            
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
        InstrumentContract _instrumentContract;

        // Constructor
        public Instrument(Pricer pricer, InstrumentContract instrumentContract)
        {
            _pricer = pricer;
            _instrumentContract = instrumentContract;
        }

        // -> METHODS 
        // --------------------
        public virtual double price(double date, int path, IScenario scenario)
        {
            return _pricer.price(date, path, _instrumentContract, scenario);
        }

        // -> ACCEORS & MODIFIERS
        // -----------------------
        public Pricer pricer { get { return _pricer; } set { _pricer = value; } }
        public InstrumentContract instrumentContract { get { return _instrumentContract; } set { _instrumentContract = value; } }
    }

    /// <summary>
    ///     Information about the contractual agreement.
    ///     Populate with NEEDED ATTRIBUTES
    /// </summary>
    
    public class InstrumentContract 
    {
        /*public enum PayoffType : uint
        {
            Call = 1,
            Put = 2,
            Forward =3,
        }*/

        string _underlyingID;   // Underlying or risk factor to use in the payoff
        double _strike;         // Strike level
        string _optionType; // Option type ("call", "put")
        double _optionDate;
        double _nominal;
        double _vol;
        double _riskFreeRate;

        public InstrumentContract(string underlyingID, string optionType, double strike, double optionDate, double nominal, double vol, double riskFreeFactor)
        {
            _underlyingID = underlyingID;
            _strike = strike;
            _optionType = optionType;
            _optionDate = optionDate;
            _nominal = nominal;
            _vol = vol;
            _riskFreeRate = riskFreeRate;
        }

        // -> ACCESORS & MODIFIERS
        // -----------------------
        public double strike { get { return _strike; } set { _strike = value; } }
        public string optionType { get { return _optionType; } set { _optionType = value; } }
        public string underlyingID { get { return _underlyingID; } set {_underlyingID = value; } }
        public double optionDate { get { return _optionDate; } set { _optionDate = value; } }
        public double nominal { get { return _nominal; } set { _nominal = value; } }
        public double vol { get { return _vol; } set { _vol = value; } }
        public double riskFreeRate { get { return _riskFreeRate; } set { _vol = value; } }


    }

    /// <summary>
    ///     Implements the logic of pricing.
    /// </summary>
    public abstract class Pricer
    {
        public abstract double price(double date, int path, InstrumentContract instrumentContract, IScenario scenario);
    }

    public class PricerCall : Pricer
    {
        public override double price(double date, int path, InstrumentContract instrumentContract, IScenario scenario)
        {   //Call Black-Scholes 
            InstrumentContract Contract = instrumentContract as InstrumentContract;
            double dcf = (Contract.optionDate - date) / 365.25;
            if (dcf < 0) return 0.0;
            double spot = scenario.getValue(Contract.underlyingID, date, path);
            double tmp = Contract.vol * Math.Sqrt(dcf);
            double d1 = (Math.Log(spot / Contract.strike) + (Contract.riskFreeRate + (Contract.vol * Contract.vol) * 0.5) * dcf) / tmp;
            double d2 = d1 - tmp;
            return (spot * Functions.N(d2)) - 
                   (Contract.strike * Math.Exp(-Contract.riskFreeRate * dcf) * Functions.N(d1));

        }
    }

    public class PricerPut : Pricer
    {
        public override double price(double date, int path, InstrumentContract instrumentContract, IScenario scenario)
        {   //Put Black-Scholes
            InstrumentContract Contract = instrumentContract as InstrumentContract;
            double dcf = (Contract.optionDate - date) / 365.25;
            if (dcf <= 0) return 0.0;
            double spot = scenario.getValue(Contract.underlyingID, date, path);
            double tmp = Contract.vol * Math.Sqrt(dcf);
            double d1 = (Math.Log(spot / Contract.strike) + (Contract.riskFreeRate + (Contract.vol * Contract.vol) * 0.5) * dcf) / tmp;
            double d2 = d1 - tmp;

            return (Contract.strike * Math.Exp(-Contract.riskFreeRate * dcf) * Functions.N(-d2)) -
                    (spot * Functions.N(-d1));
        }
    }

    public class PricerForward : Pricer
    {
        public override double price(double date, int path, InstrumentContract instrumentContract, IScenario scenario)
        {
            InstrumentContract Contract = instrumentContract as InstrumentContract;
            double dcf = (Contract.optionDate - date) / 365.25;
            if (dcf <= 0) return 0.0;
            double Kdisc = Contract.strike * Math.Exp(-Contract.riskFreeRate * dcf);
            double spot = scenario.getValue(Contract.underlyingID, date, path);
            return spot - Kdisc;
        }
    }

    // --------------------------------------------------------------------------------------




    // METRICS (HISTOGRAM, ETC ..)  -----------------------------------------------------------------------------------

    public abstract class Metric
    {
        protected static Dictionary<double, List <double> >  _NPVList;
        public abstract void addNPVToMetric(double date, int path, double NPV, IScenario scenario);
        public abstract void showMetric();

    }
    /// <summary>
    ///     Stores the Distribution for the NPV for a certain date
    /// </summary>
    public class Histogram : Metric
    {
        private int _nBars;

        public Histogram(int nBars)
        {
            _NPVList = new Dictionary<double, List <double> >();
            _nBars = nBars;
        }

        public override void addNPVToMetric(double date, int path, double NPV, IScenario scenario)
        {
            List<double> NPVdate;
            if(_NPVList.TryGetValue(date,out NPVdate) == false)
            {
                List<double> newList = new List<double>();
                newList.Add(NPV);
                _NPVList.Add(date, newList);
            }
            else
                NPVdate.Add(NPV);  
        }
        public override void showMetric()
        {
            foreach (double Date in _NPVList.Keys)
            {
                HistogramByDate(Date);
            }
        }
        public void HistogramByDate(double date)
        {
            List<double> NPVs = new List<double>(_NPVList[date]);
            double[] probability = new double[_nBars];
            NPVs.Sort();
            
            double min = NPVs[0];
            double max = NPVs[NPVs.Count - 1];
            double range = max - min;
            double step = (1.0 / (_nBars -1)) * range;
            int k = 0;
            double bound = step;
            if (max == 0 && min == 0) return;
            Console.WriteLine("Histograma NPV a la fecha {0}", date);
            Console.WriteLine("");
            
            for (int i = 0; i < NPVs.Count; i++)
            {

                while ((NPVs[i] - min) > bound)
                {
                    k++;
                    bound += step;
                }
                probability[k] += 1;
                    
            }
            bound = 0;
            for(int i = 0; i < _nBars; i++)
            {
                probability[i] /= NPVs.Count;
                Console.Write("[{0:0.00}][{1:0.00}]", bound + min,bound + step + min);
                for (int j = 0; j < Convert.ToInt16(probability[i] * 100); j++)
                    Console.Write("=");
                Console.Write(">     P = {0}", probability[i]);
                Console.Write("\n");
                bound += step;
            }
            Console.WriteLine("");
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
        List<double> _simDates;
        Scenario _scenario;
        int _nSimulations;

        // CONSTRUCTOR (TO BE IMPLEMENTED)
        // -------------------------------
        public MCEngine()
        {

            _l_Instruments = new List<Instrument>();
            _l_Metric = new List<Metric>();
            _l_Metric.Add(new Histogram(15));

            //Read Excel info and fill Attributes
            List<RiskFactor> Equities;
            double[,] CorrelationMatrix;
            _simDates = ReadSimulationDates();
            Equities = ReadEquities();
            CorrelationMatrix = ReadCorrelationMatrix();
            _l_Instruments = ReadInstruments(Equities);
            _nSimulations = ReadSimulationData();
            _engine = new Evolver(Equities, CorrelationMatrix);
            _scenario = new Scenario(Equities, _simDates[0]);
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

        public int ReadSimulationData()
        {
            List<double> Dates = new List<double>();
            ReadRange RangeHandle = new ReadRange();
            DataSet Dates_Data = RangeHandle.Read("Simulation");
            DataRowCollection Row = Dates_Data.Tables["Simulation"].Rows;
            return Convert.ToInt32(Row[0][1]);

        }

        public List<RiskFactor> ReadEquities()
        {
            List<RiskFactor> Equities = new List<RiskFactor>();
            ReadRange RangeHandle = new ReadRange();
            DataSet Equity_Data = RangeHandle.Read("EquityFeatures");
            DataRowCollection EquityRows = Equity_Data.Tables["EquityFeatures"].Rows;
            DataColumnCollection EquityColumns = Equity_Data.Tables["EquityFeatures"].Columns;
            string prueba = Convert.ToString(EquityRows[0][0]);

            foreach (DataRow DRow in EquityRows)
            {
                Equities.Add(new EquityRiskFactor(Convert.ToString(DRow[0]), Convert.ToDouble(DRow[1]), Convert.ToDouble(DRow[2]), Convert.ToDouble(DRow[3])));
            }
            return Equities;

        }

        public double[,] ReadCorrelationMatrix()
        {
            double[,] correlMatrix;
            ReadRange Matrix_Range = new ReadRange();
            DataSet Matrix_Data = Matrix_Range.Read("CorrelationMatrix");
            DataColumnCollection Column = Matrix_Data.Tables["CorrelationMatrix"].Columns;
            DataRowCollection Row = Matrix_Data.Tables["CorrelationMatrix"].Rows;

            correlMatrix = new double[Row.Count, Column.Count-1];

            for (int i = 0; i < Row.Count; i++)
            {
                for (int j = 1; j < Column.Count; j++)
                {
                   correlMatrix[i, j-1] = Convert.ToDouble(Row[i][j]);
                }
            }

            return correlMatrix;
        }

        public List<Instrument> ReadInstruments(List<RiskFactor> RFactors)
        {
            // Reading instrument contracts for our portfolio

            List<Instrument> instruments = new List<Instrument>();
            ReadRange Instrument_Range = new ReadRange();
            DataSet Instrument_Data = Instrument_Range.Read("Instruments");
            DataColumnCollection Column = Instrument_Data.Tables["Instruments"].Columns;
            DataRowCollection Rows = Instrument_Data.Tables["Instruments"].Rows;

            foreach(DataRow Row in Rows)
            {
                EquityRiskFactor RF = RFactors.Find(p => p.Equals(Convert.ToString(Row[0]))) as EquityRiskFactor;
                InstrumentContract contract = new InstrumentContract(Convert.ToString(Row[0]),
                                                                                Convert.ToString(Row[1]), 
                                                                                Convert.ToDouble(Row[2]),
                                                                                Convert.ToDouble(Row[3]), 
                                                                                Convert.ToDouble(Row[4]),
                                                                                RF.vol,
                                                                                RF.riskFreeRate);
                Pricer pricer = FactoryPricer.CreatePricer(Convert.ToString(Row[1]));
                instruments.Add(new Instrument(pricer, contract));
            }

            return instruments;
        }

        // -> METHODS
        // -------------------------------
        public void calculate()
        {
     

            // -> CREATE SCENARIO ...
            for (int j = 0; j < _nSimulations; ++j)
            {
                InitializePath(j);
                for (int i = 1; i < _simDates.Count; ++i)
                {
                    // -> EVOLVE ..
                    _engine.evolve(_simDates[i - 1], _simDates[i], j, _scenario);

                    // -> PRICE
                    double npv = 0.0;
                    for (int k = 0; k < _l_Instruments.Count; ++k)
                        npv += _l_Instruments[k].price(_simDates[i], j, _scenario);

                    // -> CALL TO METRIC TO TAKE INTO ACCOUNT THE NPV
                    foreach (Metric metric in _l_Metric)
                        metric.addNPVToMetric(_simDates[i], j, npv, _scenario);
                }
            }

        }

        public void ShowMetrics()
        {
            foreach (Metric Mtc in _l_Metric)
            {
                Mtc.showMetric();
            }
        }

        public void InitializePath(int i)
        {
            List<RiskFactor> RFList = _engine.RiskFactors;
            foreach (RiskFactor RF in RFList)
            {
                _scenario.setValue(RF.ID, _simDates[0], i, RF.spot);
            }
        }


       
    }

    public class ReadRange
    {
        public DataSet Read(string RangeName)
        {
            String sConnectionString = "Provider=Microsoft.Jet.OLEDB.4.0;Data Source="+System.Environment.CurrentDirectory+"\\Duffy_Interface.xls;Extended Properties=Excel 8.0;";
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

    public  class Functions
    {
        public static double n(double x)
        {
            double A = 1.0 / Math.Sqrt(2.0 * 3.1415);
            return A * Math.Exp(-x * x * 0.5);
        }

        public static double N(double x)
        {
            double a1 = 0.43661836;
            double a2 = -0.1201676;
            double a3 = 0.9372980;

            double k = 1.0 / (1.0 + (0.33267 * x));

            if(x >= 0.0)
            {
                return 1.0 -n(x)*(a1*k + (a2 * k * k) + (a3 *k *k *k));
            }
            else
            {
                return 1.0 -N(-x);
            }
        }

        public static void swap(ref int x,ref int y)
        {
            int aux = x;
            x = y;
            y = aux;
        }

     }

    public class FactoryPricer
    {
        public static Pricer CreatePricer(string InstrumentType)
        {
            if (InstrumentType == "Call")
                return new PricerCall();
            else if (InstrumentType == "Put")
                return new PricerPut();
            else
                return new PricerForward();
        }
    }



}
