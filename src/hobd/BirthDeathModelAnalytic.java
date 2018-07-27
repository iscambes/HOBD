package hobd;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.TreeParser;

public class BirthDeathModelAnalytic extends SpeciesTreeDistribution {

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

    public double get_c1(double b, double d, double s) {
        return Math.abs(Math.sqrt((b - d - s)*(b - d - s) + 4*b*s));
    }

    public double get_c2(double b, double d, double s, double rho, double c1) {
        return -(b - d - 2*b*rho - s)/c1;
    }

    public double get_p0(double t, double b, double d, double s, double c1, double c2) {

        return (b + d + s
                + c1*(Math.exp(-c1*t)*(1 - c2) - (1 + c2)) / (Math.exp(-c1*t)*(1 - c2) + (1 + c2)))
                / (2.0*b);
    }

    double get_q(double t, double c1, double c2) {

        return 2*(1 - c2*c2) + Math.exp(-c1*t)*Math.pow(1-c2,2.0) + Math.exp(c1*t)*Math.pow(1+c2, 2.0);
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        if (tree.getRoot().getHeight() > timeOrigin.getValue())
            return Double.NEGATIVE_INFINITY;

        double logP = 0.0;

        double c1 = get_c1(birthRate.getValue(), deathRate.getValue(), samplingRate.getValue());
        double c2 = get_c2(birthRate.getValue(), deathRate.getValue(), samplingRate.getValue(), presentSamplingProb.getValue(), c1);


        Node[] treeNodes = tree.getNodesAsArray();
        for (Node thisNode : treeNodes) {
            if (thisNode.isLeaf()) {
                if (thisNode.getHeight() > 0.0 || presentSamplingProb.getValue() == 0.0) {
                    // psi-sample

                    logP = logP + Math.log(samplingRate.getValue())
                            + Math.log(get_q(thisNode.getHeight(), c1, c2));
//                            + Math.log(get_p0(thisNode.getHeight(),
//                                    birthRate.getValue(),
//                                    deathRate.getValue(),
//                                    samplingRate.getValue(), c1, c2));
                } else {
                    // sample at present day

                    logP = logP + Math.log(4.0 * presentSamplingProb.getValue());
                }
            } else {
                // coalescence node

                logP = logP + Math.log(birthRate.getValue()) - Math.log(get_q(thisNode.getHeight(), c1, c2));
            }
        }

        logP = logP - Math.log(get_q(timeOrigin.getValue(), c1, c2));

        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }


    public static void main(String[] args) {

        double b = 1.0, d = 0.5, s = 0.6, rho = 0.2, t0= 5.0;

        BirthDeathModelAnalytic bdmodel = new BirthDeathModelAnalytic();

        TreeParser tree = new TreeParser("(A:1.0, B:1.0):0.0;");
        bdmodel.initByName("birthRate", new RealParameter(String.valueOf(b)),
                "deathRate", new RealParameter(String.valueOf(d)),
                "samplingRate", new RealParameter(String.valueOf(s)),
                "presentSamplingProb", new RealParameter(String.valueOf(rho)),
                "timeOrigin", new RealParameter(String.valueOf(t0)),
                "tree", tree);


        double c1 = bdmodel.get_c1(b, d, s);
        double c2 = bdmodel.get_c2(b, d, s, rho, c1);

        double res = bdmodel.get_q(5.0, c1, c2);
        double tete = bdmodel.get_p0(5.0, b, d, s, c1, c2);


        System.out.println("q = " + res);
    }
}
