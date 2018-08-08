package hobd.simulator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class HOBDTrajectory {
    protected List<HOBDEvent> eventList;
    protected List<HOBDState> stateList;

    protected List<HOBDEvent> samplingEvents;

    protected boolean stateListDirty = true;

    public HOBDTrajectory() {
        stateListDirty = false;
    }

    public HOBDTrajectory(List<HOBDEvent> eventList) {
        this.eventList = new ArrayList<>(eventList);
        stateListDirty = true;
    }

    public void setEventList(List<HOBDEvent> eventList) {
        this.eventList = eventList;
        stateListDirty = true;
    }

    public List<HOBDEvent> getEventList() {
        return eventList;
    }

    public List<HOBDState> getStateList() {
        updateStateList();

        return stateList;
    }

    private void updateStateList() {
        if (!stateListDirty)
            return;

        stateListDirty = false;

        if (eventList == null) {
            stateList = null;
            return;
        }

        stateList = new ArrayList<>();
        samplingEvents = new ArrayList<>();

        long currentPopSize = 1;

        stateList.add(new HOBDState(0.0, currentPopSize));

        for (HOBDEvent event : eventList) {
            currentPopSize += event.getDelta();
            stateList.add(new HOBDState(event.time, currentPopSize));

            if (event.isSample() && event.count>0)
                samplingEvents.add(event);
        }
    }

    public int getSampleCount() {
        int count = 0;

        for (HOBDEvent event : getSamplingEvents())
            count += event.count;

        return count;
    }

    public List<HOBDEvent> getSamplingEvents() {
        updateStateList();

        return samplingEvents;
    }

    public void print(PrintStream ps) {

        updateStateList();

        ps.println("time\tevent\tpopSize");

        for (int i=0; i<eventList.size(); i++)
            ps.println(eventList.get(i).time + "\t"
                    + eventList.get(i).type + "\t"
                    + stateList.get(i).popSize);
    }
}
