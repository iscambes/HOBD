package hobd.simulator;

public class HOBDEvent {
    public enum HOBDEventType {BIRTH, BURST, DEATH, PSISAMPLE, RHOSAMPLE};
    public HOBDEventType type;
    public double time;
    public int count;

    public long getDelta() {
        switch(type) {
            case BIRTH:
                return 1;
            case BURST:
                return count;
            case DEATH:
            case PSISAMPLE:
            case RHOSAMPLE:
                return -count;
            default:
                throw new IllegalArgumentException("Unknown event type.");
        }
    }

    public boolean isSample() {
        return type == HOBDEventType.PSISAMPLE || type == HOBDEventType.RHOSAMPLE;
    }

    @Override
    public String toString() {
        return type + " (count " + count + ") at time " + time;
    }
}
