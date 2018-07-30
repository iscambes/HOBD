package hobd;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.TreeParser;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaIntegrator;

public class BirthDeathModelNumeric extends SpeciesTreeDistribution {

    public Input<RealParameter> birthRateInput = new Input<>("birthRate",
            "Per-individual birth rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> deathRateInput = new Input<>("deathRate",
            "Per-individual death rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingRateInput = new Input<>("samplingRate",
            "Per-individual sampling rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> presentSamplingProbInput = new Input<>("presentSamplingProb",
            "Per-individual present-day sampling probability.", Input.Validate.REQUIRED);

    public Input<RealParameter> timeOriginInput = new Input<>("timeOrigin",
            "Per-individual time Origin.", Input.Validate.REQUIRED);



    RealParameter birthRate, deathRate, samplingRate, presentSamplingProb, timeOrigin;

    @Override
    public void initAndValidate() {
        birthRate = birthRateInput.get();
        deathRate = deathRateInput.get();
        samplingRate = samplingRateInput.get();
        presentSamplingProb = presentSamplingProbInput.get();
        timeOrigin = timeOriginInput.get();


    }

    class ODEp0 implements FirstOrderDifferentialEquations {

        double lambda, mu, psi, rho;

        public ODEp0(double lambda, double mu, double psi, double rho) {
            this.lambda = lambda;
            this.mu = mu;
            this.psi = psi;
            this.rho = rho;
        }

        @Override
        public int getDimension() {
            return 1;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot)
                throws MaxCountExceededException, DimensionMismatchException {

            double p0 = y[0];

            yDot[0] = -(lambda + mu + psi)*p0 + mu + lambda*p0*p0;
        }
    }

    public double get_p0(double t) {

        ODEp0 p0equations = new ODEp0(birthRate.getValue(), deathRate.getValue(),
                samplingRate.getValue(), presentSamplingProb.getValue());

        // AbstractIntegrator integrator = new EulerIntegrator(0.001);
        AbstractIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.001);

        double [] state = {1 - presentSamplingProb.getValue()};

        if (t>0) {
            integrator.integrate(p0equations, 0, state, t, state);
        }

        return state[0];
    }

    public class ODEgep0 implements FirstOrderDifferentialEquations {

        double lambda, mu, psi, rho;

        public ODEgep0(double lambda, double mu, double psi, double rho) {
            this.lambda = lambda;
            this.mu = mu;
            this.psi = psi;
            this.rho = rho;
        }


        @Override
        public int getDimension() {
            return 2;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {

            double p0 = y[0];
            double ge = y[1];

            double p0Dot = ;
            double geDot = ;

            yDot[0] = p0Dot;
            yDot[1] = geDot;
        }
    }

    double get_ge(double t, Node node) {


        return 0;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        double logP = 0.0;

        logP = get_ge(timeOrigin.getValue(), tree.getRoot());

        return logP;
    }

    public static void main(String[] args) {

        double d = 0.5, s = 0.6, rho = 0.2, t0 = 5.0, t = 2.0;

        TreeParser tree = new TreeParser("(A:1.0, B:1.0):0.0;");

        for (double b = 0.5; b<= 1.5; b += 0.1) {

            BirthDeathModelNumeric bdmodelNumeric = new BirthDeathModelNumeric();
            bdmodelNumeric.initByName("birthRate", new RealParameter(String.valueOf(b)),
                    "deathRate", new RealParameter(String.valueOf(d)),
                    "samplingRate", new RealParameter(String.valueOf(s)),
                    "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                    "timeOrigin", new RealParameter(String.valueOf(t0)),
                    "tree", tree);

            double p0numeric = bdmodelNumeric.get_p0(t);

            BirthDeathModelAnalytic bdmodelAnalytic = new BirthDeathModelAnalytic();
            bdmodelAnalytic.initByName("birthRate", new RealParameter(String.valueOf(b)),
                    "deathRate", new RealParameter(String.valueOf(d)),
                    "samplingRate", new RealParameter(String.valueOf(s)),
                    "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                    "timeOrigin", new RealParameter(String.valueOf(t0)),
                    "tree", tree);

            double c1 = bdmodelAnalytic.get_c1(b, d, s);
            double c2 = bdmodelAnalytic.get_c2(b, d, s, rho, c1);

            double p0analytic = bdmodelAnalytic.get_p0(t, b, d, s, c1, c2);

            System.out.println(b + "\t" + p0analytic + "\t" + p0numeric);
        }
    }
}
