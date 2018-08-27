package hobd;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.math.GammaFunction;
import beast.math.distributions.Gamma;
import beast.util.TreeParser;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import pitchfork.Pitchforks;

import java.util.ArrayList;
import java.util.List;



public class HOBDModel extends SpeciesTreeDistribution {

    public Input<HOBDModelParams> modelParamsInput = new Input<>("modelParams",
            "Higher-order birth-death model parameters.",
            Input.Validate.REQUIRED);

    HOBDModelParams modelParams;

    @Override
    public void initAndValidate() {
        modelParams = modelParamsInput.get();
    }

    class ODEp0 implements FirstOrderDifferentialEquations {

        double lambda, mu, psi, rho, meanBurstSize, burstRate;

        public ODEp0(double lambda, double mu, double psi, double rho, double meanBurstSize, double burstRate) {
            this.lambda = lambda;
            this.mu = mu;
            this.psi = psi;
            this.rho = rho;
            this.meanBurstSize = meanBurstSize;
            this.burstRate = burstRate;


        }

        @Override
        public int getDimension() {
            return 1;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot)
                throws MaxCountExceededException, DimensionMismatchException {

            double p0 = y[0];

            yDot[0] = -(lambda + mu + psi + burstRate)*p0 + mu + (lambda + burstRate * Math.exp((meanBurstSize-1)*(p0-1))*p0*p0);


        }
    }


    public double get_p0(double t) {

        ODEp0 p0equations = new ODEp0(modelParams.getBirthRate(), modelParams.getDeathRate(),
                modelParams.getSamplingRate(), modelParams.getPresentSamplingProb(),
                modelParams.getBurstRate(), modelParams.getMeanBurstSize());

        //AbstractIntegrator integrator = new EulerIntegrator(0.02);
        AbstractIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.001);

        double [] state = {1 - modelParams.getPresentSamplingProb()};

        if (t>0) {
            integrator.integrate(p0equations, 0, state, t, state);
        }

        return state[0];
    }


    public class ODEgep0 implements FirstOrderDifferentialEquations {

        double lambda, mu, psi, rho, meanBurstSize, burstRate;

        public ODEgep0(double burstRate, double lambda, double mu, double psi, double rho, double meanBurstSize) {
            this.lambda = lambda;
            this.mu = mu;
            this.psi = psi;
            this.rho = rho;
            this.meanBurstSize = meanBurstSize;
            this.burstRate = burstRate;
        }
        @Override
        public int getDimension() { return 2;
        }

        @Override
        public void computeDerivatives(double t, double[] y, double[] yDot) {

            double p0 = y[0];
            double ge = y[1];

            ////ECUACIONES DIFERENCIALES.

            double p0Dot = -(lambda + mu + psi + burstRate)*p0 + mu + (lambda + burstRate * Math.exp((meanBurstSize-1)*(p0-1)))*p0*p0;
            double geDot = -(lambda + mu + psi + burstRate)*ge + (2*lambda*p0*ge)
                    + burstRate*Math.exp((meanBurstSize-1)*(p0-1))*(2 + (meanBurstSize-1)*p0)*p0*ge;



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
                state[1] = modelParams.getSamplingRate();

            } else {

                // rho-sampling

                state[0] = 1-modelParams.getPresentSamplingProb();
                state[1] = modelParams.getPresentSamplingProb();
            }

        } else {

            //internal node

            List<Node> children = Pitchforks.getLogicalChildren(node);
            List<double[]> childStates = new ArrayList<>();

            int k = children.size();

            for (Node child : children) {
                childStates.add(get_gep0(node.getHeight(), child));
            }

            double term1;
            if (k==2) {
                term1 = modelParams.getBirthRate() * childStates.get(0)[1] * childStates.get(1)[1];
            } else {
                term1 = 0.0;
            }

            double p0 = childStates.get(0)[0];

            double acc = Math.exp((modelParams.getMeanBurstSize()-1)*p0);

            for (int n=1; n<=k-2; n++) {
                acc -= Math.exp((n-1)*(Math.log(modelParams.getMeanBurstSize()-1) + Math.log(p0))
                        - GammaFunction.lnGamma(n-1 + 1));
            }

            double term2 = modelParams.getBurstRate()*Math.exp(-(modelParams.getMeanBurstSize()-1))
                    * 1.0/Math.pow(p0, k-2)
                    * acc;

            for (double[] childState : childStates)
                term2 *= childState[1];

            state[0] = p0;
            state[1] = term1 + term2;
        }


        ODEgep0 gep0equations = new ODEgep0(modelParams.getBirthRate(),modelParams.getDeathRate(),
                modelParams.getSamplingRate(), modelParams.getPresentSamplingProb(),
                modelParams.getBurstRate(), modelParams.getMeanBurstSize());

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

        logP = Math.log(get_ge(modelParams.getOriginTime(), tree.getRoot())); // remember that    //f[T |λ,μ,ψ, tor = x0] = ge(tor)//

        return logP;
    }


    public static void main(String[] args) {

        double b=1.0, d = 0.5, s = 0.6, rho = 0.2, t0 = 5.0, m = 20, j = 10;

        TreeParser tree = new TreeParser("(A:1.0,B:1.0):1.0;");   // con sto del tree pasrser podiamos ir analizando como de uno en un edge. la parte de arriba es para varios (dos creo) y lo de abajo que esta con // para UNO SOLO

//        TreeParser tree = new TreeParser("A:1.0;");


        //SUPERSPREADING MODEL

        HOBDModel HOBDModel = new HOBDModel();
        HOBDModel.initByName("birthRate", new RealParameter(String.valueOf(b)),
                "deathRate", new RealParameter(String.valueOf(d)),
                "samplingRate", new RealParameter(String.valueOf(s)),
                "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                "timeOrigin", new RealParameter(String.valueOf(t0)),
                "mean", new RealParameter(String.valueOf(m)),
                "superRate", new RealParameter(String.valueOf(j)),
                "tree", tree);

        double p0HOBDModel = HOBDModel.get_p0(t0);
        double geHOBDMODEL = HOBDModel.get_ge(t0, tree.getRoot());

        //NUMERICALMODEL

        BirthDeathModelNumeric bdmodelNumeric = new BirthDeathModelNumeric();
        bdmodelNumeric.initByName("birthRate", new RealParameter(String.valueOf(b)),
                "deathRate", new RealParameter(String.valueOf(d)),
                "samplingRate", new RealParameter(String.valueOf(s)),
                "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                "timeOrigin", new RealParameter(String.valueOf(t0)),
                "tree", tree);

        double p0numeric = bdmodelNumeric.get_p0(t0);
        double geNumeric = bdmodelNumeric.get_ge(t0, tree.getRoot());

        //ANALYTIC MODEL

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



        System.out.println("p0: " + p0analytic + "\t" + p0HOBDModel + "\t" + p0numeric);
        System.out.println("ge: " + geAnalytic + "\t" + geHOBDMODEL + "\t" + geNumeric);






    }
}