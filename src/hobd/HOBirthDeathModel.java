package hobd;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.TreeParser;

public class HOBirthDeathModel extends SpeciesTreeDistribution {

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

    double c1() {
        return Math.abs(Math.sqrt(Math.pow(birthRate.getValue()-deathRate.getValue()-samplingRate.getValue(), 2.0)
                + 4*birthRate.getValue()*samplingRate.getValue()));
    }

    double c2() {
        return Math.abs(-((birthRate.getValue()-deathRate.getValue()-(2*birthRate.getValue()*presentSamplingProb.getValue())-samplingRate.getValue())/c1()));

    }

    double p0(double t) {
        return Math.abs(((birthRate.getValue()+deathRate.getValue()+samplingRate.getValue())+c1()*((-(1+c2())+(1-c2())*Math.exp(-c1()*t))
                /((1+c2())+(1-c2())*Math.exp(-c1()*t))))/(2*birthRate.getValue()));
    }

    double q(double t) {
        return Math.abs((2*(1-(Math.pow(c2(), 2.0))))+((Math.exp(-c1()*t))*(Math.pow(1-c2(), 2.0)))+(((Math.exp(c1()*t))*(Math.pow(1+c2(), 2.0)))));
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        logP = 0.0;

        Node[] treeNodes = tree.getNodesAsArray();
        for (Node thisNode : treeNodes) {
            if (thisNode.isLeaf()) {
                if (thisNode.getHeight() > 0.0) {
                    // psi-sample

                    logP = logP + Math.log(samplingRate.getValue()) + Math.log(p0(thisNode.getHeight()) * q(thisNode.getHeight()));
                } else {
                    // sample at present day

                    logP = logP + Math.log((4.0 * presentSamplingProb.getValue()));
                }
            } else {
                // coalescence node

                logP = logP + Math.log(birthRate.getValue()) + Math.log(1.0 / q(thisNode.getHeight()));
            }
        }

        logP = logP - Math.log(q(timeOrigin.getValue()));

        return logP;
    }

    public static void main(String[] args) {

        HOBirthDeathModel hobd = new HOBirthDeathModel();

        TreeParser tree = new TreeParser("(A:1.0, B:1.0):0.0;");
        hobd.initByName("birthRate", new RealParameter("1.0"),
                "deathRate", new RealParameter("0.5"),
                "samplingRate", new RealParameter("0.6"),
                "presentSamplingProb", new RealParameter("0.2"),
                "timeOrigin", new RealParameter("5.0"),
                "tree", tree);

        System.out.println("c1 = " + hobd.c1());
        System.out.println("c2 = " + hobd.c2());
        System.out.println("p0(t) = " + hobd.p0(1.0));
        System.out.println("q(t) = " + hobd.q(1.0));
    }
}
