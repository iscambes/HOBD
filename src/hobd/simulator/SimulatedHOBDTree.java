package hobd.simulator;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.util.Randomizer;
import hobd.HOBDModelParams;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class SimulatedHOBDTree extends Tree {

    public Input<HOBDModelParams> paramsInput = new Input<>("modelParams",
            "HOBD model params.", Input.Validate.REQUIRED);

    public Input<Integer> minSampleCountInput = new Input<>("minSampleCount",
            "Minimum sample count.", 1);

    HOBDModelParams params;
    int minSampleCount;

    @Override
    public void initAndValidate() {

        params = paramsInput.get();
        minSampleCount = minSampleCountInput.get();

        super.initAndValidate();

        simulate();
    }

    void simulate() {

        HOBDTrajectory traj;
        do {
            traj = new SimulatedHOBDTrajectory(params);
        } while (traj.getSampleCount() < minSampleCount);

        List<HOBDEvent> eventList = traj.getEventList();
        List<HOBDState> stateList = traj.getStateList();

        Collections.reverse(eventList);
        Collections.reverse(stateList);

        List<Node> activeLineages = new ArrayList<>();

        int samplesRemaining = traj.getSampleCount();
        int nextInternalNodeNr = traj.getSampleCount();

        int eventIdx = 0;
        while (samplesRemaining>0 || activeLineages.size()>1) {

            HOBDEvent event = eventList.get(eventIdx);
            HOBDState state = stateList.get(eventIdx);

            int n = activeLineages.size();
            double I = state.popSize;
            double age = params.getOriginTime() - event.time;

            double pCoal;
            Node parent;

            switch(event.type) {
                case BIRTH:
                    pCoal = n*(n-1.0)/(I*(I-1.0));
                    if (Randomizer.nextDouble()>pCoal)
                        break;

                    Node left = activeLineages.get(Randomizer.nextInt(n));
                    activeLineages.remove(left);
                    Node right = activeLineages.get(Randomizer.nextInt(n-1));
                    activeLineages.remove(right);

                    parent = new Node();
                    parent.setNr(nextInternalNodeNr++);
                    parent.addChild(left);
                    parent.addChild(right);
                    parent.setHeight(age);

                    activeLineages.add(parent);
                    break;

                case BURST:
                    if (n<2)
                        break;

                    double[] p = new double[n-1];

                    int burstSize = event.count;

                    pCoal = 0.0;
                    for (int k=2; k<=n; k++) {
                        p[k-2] = Math.exp(Binomial.logChoose(burstSize, k) + Binomial.logChoose(n, k) - Binomial.logChoose((int)Math.round(I), k));
                        pCoal += p[k-2];
                    }

                    if (Randomizer.nextDouble()>pCoal)
                        break;

                    double u = Randomizer.nextDouble()*pCoal;
                    int k;
                    for (k=2; k<=n; k++) {
                        if (u<p[k-2])
                            break;

                        u -= p[k-2];
                    }

                    System.out.println("Burst size " + burstSize + " k=" + k);

                    parent = new Node();
                    parent.setNr(nextInternalNodeNr++);
                    parent.setHeight(age);

                    for (int i=0; i<k; i++) {
                        Node child = activeLineages.get(Randomizer.nextInt(n));
                        activeLineages.remove(child);
                        n -= 1;

                        parent.addChild(child);
                    }

                    activeLineages.add(parent);
                    break;

                case PSISAMPLE:
                case RHOSAMPLE:
                    for (int i=0; i<event.count; i++) {
                        Node sampleNode = new Node();
                        sampleNode.setNr(samplesRemaining-1);
                        sampleNode.setHeight(age);
                        activeLineages.add(sampleNode);

                        samplesRemaining -= 1;
                    }
                    break;

                default:
                    break;
            }

            eventIdx += 1;
        }

        assignFromWithoutID(new Tree(activeLineages.get(0)));
    }

    public static void main(String[] args) {

//        for (int s=1; s<100; s++) {
//            Randomizer.setSeed(s);
//
//            System.out.println("** " + s + " **");
        Randomizer.setSeed(3);

            SimulatedHOBDTree tree = new SimulatedHOBDTree();

            HOBDModelParams params = new HOBDModelParams(1.3, 0.8,
                    0.1, 10, 0.3, 0.1,
                    3.0);

            tree.initByName("modelParams", params);

            System.out.println(tree);
//        }
    }
}
