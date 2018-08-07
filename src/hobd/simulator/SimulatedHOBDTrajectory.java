package hobd.simulator;

import beast.util.Randomizer;
import hobd.HOBDModelParams;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class SimulatedHOBDTrajectory extends HOBDTrajectory {

    public HOBDModelParams params;

    public SimulatedHOBDTrajectory(HOBDModelParams params) {

        this.params = params;

        simulate();
    }

    public void simulate() {

        eventList = new ArrayList<>();

        HOBDState currentState = new HOBDState(0.0, 1);

        Map<HOBDEvent.HOBDEventType,Double> propensities = new HashMap<>();

        while (true) {

            currentState.updatePropensities(params, propensities);

            double totalProp = 0.0;
            for (double val : propensities.values())
                totalProp += val;

            if (totalProp>0.0)
                currentState.time += Randomizer.nextExponential(totalProp);
            else
                currentState.time = Double.POSITIVE_INFINITY;

            if (currentState.time > params.getOriginTime()){
                break;
            }

            HOBDEvent event = new HOBDEvent();
            event.time = currentState.time;

            double u = Randomizer.nextDouble()*totalProp;

            for (HOBDEvent.HOBDEventType thisEventType : propensities.keySet()) {
                if (u < propensities.get(thisEventType)) {
                    event.type = thisEventType;
                    break;
                } else {
                    u -= propensities.get(thisEventType);
                }
            }

            if (event.type == null)
                throw new IllegalStateException("Reaction selection loop fell through.");

            if (event.type == HOBDEvent.HOBDEventType.BURST)
                event.count = 2+(int)Randomizer.nextPoisson(params.getBurstSize());
            else
                event.count = 1;

            eventList.add(event);
            currentState.popSize += event.getDelta();
        }

        // Rho sampling

        HOBDEvent rhoSamplingEvent = new HOBDEvent();
        rhoSamplingEvent.time = params.getOriginTime();
        rhoSamplingEvent.type = HOBDEvent.HOBDEventType.RHOSAMPLE;
        rhoSamplingEvent.count = 0;
        for (long i=0; i<currentState.popSize; i++) {
            if (Randomizer.nextDouble() < params.getPresentSamplingProb())
                rhoSamplingEvent.count += 1;
        }

        eventList.add(rhoSamplingEvent);

        // Final marker event

        HOBDEvent finalEvent = new HOBDEvent();
        finalEvent.time = params.getOriginTime();
        finalEvent.type = HOBDEvent.HOBDEventType.END;
        finalEvent.count = 1;

        eventList.add(finalEvent);

        stateListDirty = true;
    }

    public static void main(String[] args) {

        HOBDModelParams params = new HOBDModelParams(1.2, 1.0,
                0.2, 10, 0.1, 0.1,
                3.0);

        SimulatedHOBDTrajectory traj = new SimulatedHOBDTrajectory(params);

        traj.print(System.out);
    }
}
