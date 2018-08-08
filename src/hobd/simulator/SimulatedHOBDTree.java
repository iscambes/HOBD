package hobd.simulator;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
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

            int k = activeLineages.size();
            double I = state.popSize;
            double age = params.getOriginTime() - event.time;

            switch(event.type) {
                case BIRTH:
                    double pCoal = k*(k-1.0)/(I*(I-1.0));
                    if (Randomizer.nextDouble()>pCoal)
                        break;

                    Node left = activeLineages.get(Randomizer.nextInt(k));
                    activeLineages.remove(left);
                    Node right = activeLineages.get(Randomizer.nextInt(k-1));
                    activeLineages.remove(right);

                    Node parent = new Node();
                    parent.setNr(nextInternalNodeNr++);
                    parent.addChild(left);
                    parent.addChild(right);
                    parent.setHeight(age);

                    activeLineages.add(parent);
                    break;

                case BURST:

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

        Randomizer.setSeed(15);

        SimulatedHOBDTree tree = new SimulatedHOBDTree();

        HOBDModelParams params = new HOBDModelParams(1.3, 0.8,
                0.0, 10, 0.3, 0.1,
                3.0);

        tree.initByName("modelParams", params);

        System.out.println(tree);
    }
}
