package hobd.simulator;

import hobd.HOBDModelParams;

import java.util.HashMap;
import java.util.Map;

public class HOBDState {
    public double time;
    public double popSize;

    public HOBDState(double time, long popSize) {
        this.time = time;
        this.popSize = popSize;
    }

    public void updatePropensities(
            HOBDModelParams params, Map<HOBDEvent.HOBDEventType, Double> propensities) {

        propensities.clear();

        propensities.put(HOBDEvent.HOBDEventType.BIRTH,
                params.getBirthRate()*popSize);

        propensities.put(HOBDEvent.HOBDEventType.BURST,
                params.getBurstRate()*popSize);

        propensities.put(HOBDEvent.HOBDEventType.DEATH,
                params.getDeathRate()*popSize);

        propensities.put(HOBDEvent.HOBDEventType.PSISAMPLE,
                params.getSamplingRate()*popSize);
    }

    @Override
    public String toString() {
        return "PopSize: " + popSize + " Time: " + time;
    }
}
