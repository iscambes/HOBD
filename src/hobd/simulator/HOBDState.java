package hobd.simulator;

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
            HOBDParams params, Map<HOBDEvent.HOBDEventType, Double> propensities) {

        propensities.clear();

        propensities.put(HOBDEvent.HOBDEventType.BIRTH,
                params.birthRate*popSize);

        propensities.put(HOBDEvent.HOBDEventType.BURST,
                params.burstRate*popSize);

        propensities.put(HOBDEvent.HOBDEventType.DEATH,
                params.deathRate*popSize);

        propensities.put(HOBDEvent.HOBDEventType.PSISAMPLE,
                params.psiSampleRate*popSize);
    }

    public Map<HOBDEvent.HOBDEventType, Double> getPropensities(HOBDParams params) {
        Map<HOBDEvent.HOBDEventType, Double> propensities = new HashMap<>();
        updatePropensities(params, propensities);

        return propensities;
    }
}
