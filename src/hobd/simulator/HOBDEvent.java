package hobd.simulator;

public class HOBDEvent {
    public enum HOBDEventType {BIRTH, BURST, DEATH, PSISAMPLE, RHOSAMPLE, END};
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
            case END:
                return 0;
            default:
                throw new IllegalArgumentException("Unknown event type.");
        }
    }

    public boolean isSample() {
        return type == HOBDEventType.PSISAMPLE || type == HOBDEventType.RHOSAMPLE;
    }
}
