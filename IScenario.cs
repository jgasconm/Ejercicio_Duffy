using System;
using System.Collections.Generic;
using System.Text;
// A VER SI SALE ESTE COMENTARIO TONTO DE ANTONIO
namespace DuffyExercise
{


    // BROWNINA MANAGER ---------------------------------------------------------------------------
    /// <summary>
    ///     Brownian Dispatcher. Given certain logic it delivers brownian motion to every risk Factor.
    ///     Implement the dispatch method. klajasdfadsfasdfafa
    /// </summary>
    public class BrownianManager
    {
        /// <summary>
        ///     Implements the logic for dispatching.
        /// </summary>
        public void initialize(List<RiskFactor> riskFactors)
        {
        }

        public double[] dispatch(string sID, double[] brownianVector)
        {
            throw(new ApplicationException("Method not implemented"));
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

    }

    /// <summary>
    ///     Class implementing the Scenario Interface
    /// </summary>
    public class Scenario : IScenario
    {

        #region IScenario Members

        public double getValue(string ID, double date, int path)
        {
            throw new Exception("The method or operation is not implemented.");
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
        protected RiskFactor()
        {
        }

        // -> ACCESORS
        // ------------------
        public string ID { get { return _ID; } }


        // -> ABSTRACT METHODS
        //---------------------------------------
        public void evolve(double dateFrom, double dateTo, int path, IScenario scenario, 
                                    double[] correlatedBrownians, BrownianManager brownianDispacther)
        {
            evolve(dateFrom, dateTo, path, scenario, brownianDispacther.dispatch(ID, correlatedBrownians));
        }
        public abstract void evolve(double dateFrom, double dateTo, int path,
                                    IScenario scenario, double[] correlatedBrownians);
        /// <summary>
        ///     Returns the number of Brownian Motions needed to evolve the Risk Factor
        /// </summary>
        public abstract int size(double date);

    }


    /// <summary>
    ///     Equity Dynamics: Implement Contructor, Abstract Method ...
    /// </summary>
    public class EquityRiskFactor : RiskFactor
    {
        // -> METHODS
        // --------------------
        public override void evolve(double dateFrom, double dateTo, int path, 
                                    IScenario scenario, double[] correlatedBrownians)
        {
            throw new Exception("The method or operation is not implemented.");
        }

        public override int size(double date)
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

    // --------------------------------------------------------------------------------------



    // MONTE CARLO ENGINE -------------------------------------------------------------------
    public abstract class Correlator
    {
        public abstract double[] correlate(double dateFrom, double dateTo, double[] independentBM);
    }


    /// <summary>
    ///     Class that evolves every Risk Factor .... 
    /// </summary>
    public class Evolver
    {
        // -> ATTRIBUTES
        // ---------------------
        List<RiskFactor>                        _l_RiskFactors;
        IGaussianGenerator                      _gaussianGenerator;
        BrownianManager                         _brownianManager;
        Correlator                              _correlator;


        // -> METHODS
        // ---------------------
        
        /// <summary>
        ///     Implement: evolve every Risk Factor from dateFrom to dateTo.
        /// </summary>
        public void evolve(double dateFrom, double dateTo, int path, IScenario scenario)
        {



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
        Pricer                          _pricer;
        InstrumentContract              _InstrumentContract;

        // -> METHODS 
        // --------------------
        public virtual double price(double date, int path, IScenario scenario)
        {
            return _pricer.price(date, path, _InstrumentContract, scenario);
        }


        // -> ACCEORS & MODIFIERS
        // -----------------------
        public Pricer pricer { get { return _pricer; } set { _pricer = value; } }
        public InstrumentContract instrumentContract{get{return _InstrumentContract;}set{_InstrumentContract = value;}}
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
        Evolver                                             _engine;
        List<Instrument>                                    _l_Instruments;
        List<Metric>                                        _l_Metric;

        // CONSTRUCTOR (TO BE IMPLEMENTED)
        // -------------------------------




        // -> METHODS
        // -------------------------------
        public void calculate(List<double> dates, int noSim)
        {
            Scenario scenario =  new Scenario();

            // -> CREATE SCENARIO ...
            for (int j = 0; j < noSim; ++j)
            {
                for (int i = 0; i < dates.Count; ++i)
                {
                    // -> EVOLVE ..
                    _engine.evolve(dates[i - 1], dates[i], j, scenario);

                    // -> PRICE
                    double npv = 0.0;
                    for (int k = 0; k < _l_Instruments.Count; ++k)
                        npv += _l_Instruments[k].price(dates[i], j,  scenario);

                    // -> CALL TO METRIC TO TAKE INTO ACCOUNT THE NPV
                    foreach(Metric metric in _l_Metric)
                        metric.addNPVToMetric(dates[i], j, npv, scenario);
                }
            }

        }       
    }

    // --------------------------------------------------------------------------------------


}
