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

            double p0Dot = -(lambda + mu + psi)*p0 + mu + lambda*p0*p0;
            double geDot = -(lambda + mu + psi)*ge + (2*lambda*p0*ge);


            yDot[0] = p0Dot;
            yDot[1] = geDot;
        }
    }

    double[] get_gep0(double t, Node node) {
        // here we have to implement recursive formula


        double[] state = new double[2];


        if (node.isLeaf()) {

            // leaf node

            if (node.getHeight() > 0.0) {

                // psi-sampling

                state[0] = get_p0(node.getHeight());
                state[1] = samplingRate.getValue();

            } else {
                // rho-sampling

                state[0] = 1-presentSamplingProb.getValue();
                state[1] = presentSamplingProb.getValue();

            }
        } else {
            //internal node

            double[] stateLeft = get_gep0(node.getHeight(), node.getLeft());
            double[] stateRight = get_gep0(node.getHeight(), node.getRight());

            state[0] = stateLeft[0]; // p0 at start of edge. we use just one of the p0 , cause are the same!!
            state[1] = birthRate.getValue() * stateLeft[1] * stateRight[1];
        }



        ODEgep0 gep0equations = new ODEgep0(birthRate.getValue(), deathRate.getValue(),
                samplingRate.getValue(), presentSamplingProb.getValue());

        // AbstractIntegrator integrator = new EulerIntegrator(0.001);
        AbstractIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.001);

        // Integrate ge and p0 along edge
        integrator.integrate(gep0equations, node.getHeight(), state, t, state);

        return state; // this is going to give us the state of one node (if it is a extant or extinct leaf, a internal node...)
    }
        // once we know the state of a node, let's calculate the ge.

    double get_ge(double t, Node node) {
        return get_gep0(t, node)[1];
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        double logP = 0.0;

        logP = Math.log(get_ge(timeOrigin.getValue(), tree.getRoot())); // remember that    //f[T |λ,μ,ψ, tor = x0] = ge(tor)//

        return logP;
    }

    public static void main(String[] args) {

        double b=1.0, d = 0.5, s = 0.6, rho = 0.2, t0 = 5.0, t = 2.0;

        TreeParser tree = new TreeParser("(A:1.0,B:1.0):1.0;");
//        TreeParser tree = new TreeParser("A:1.0;");

        BirthDeathModelNumeric bdmodelNumeric = new BirthDeathModelNumeric();
        bdmodelNumeric.initByName("birthRate", new RealParameter(String.valueOf(b)),
                "deathRate", new RealParameter(String.valueOf(d)),
                "samplingRate", new RealParameter(String.valueOf(s)),
                "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                "timeOrigin", new RealParameter(String.valueOf(t0)),
                "tree", tree);

        double p0numeric = bdmodelNumeric.get_p0(t0);
        double geNumeric = bdmodelNumeric.get_ge(t0, tree.getRoot());

        BirthDeathModelAnalytic bdmodelAnalytic = new BirthDeathModelAnalytic();
        bdmodelAnalytic.initByName("birthRate", new RealParameter(String.valueOf(b)),
                "deathRate", new RealParameter(String.valueOf(d)),
                "samplingRate", new RealParameter(String.valueOf(s)),
                "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                "timeOrigin", new RealParameter(String.valueOf(t0)),
                "tree", tree);

        double c1 = bdmodelAnalytic.get_c1(b, d, s);
        double c2 = bdmodelAnalytic.get_c2(b, d, s, rho, c1);

        double p0analytic = bdmodelAnalytic.get_p0(t0, b, d, s, c1, c2);
        double geAnalytic = Math.exp(bdmodelAnalytic.calculateTreeLogLikelihood(tree));

        System.out.println("p0: " + p0analytic + "\t" + p0numeric);
        System.out.println("ge: " + geAnalytic + "\t" + geNumeric);

    }
}
