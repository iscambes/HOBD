package hobd;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public class BirthDeathModelParams extends CalculationNode {

    public Input<RealParameter> birthRateInput = new Input<>("birthRate",
            "Per-individual birth rate", Input.Validate.REQUIRED);

    public Input<RealParameter> deathRateInput = new Input<>("deathRate",
            "Per-individual death rate", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingRateInput = new Input<>("samplingRate",
            "Per-individual sampling rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> presentSamplingProbInput = new Input<>("presentSamplingProb",
            "Per-individual present-day sampling prob.", Input.Validate.REQUIRED);

    public Input<RealParameter> originTimeInput = new Input<>("originTime",
            "Time between present and origin of process.",
            Input.Validate.REQUIRED);

    public BirthDeathModelParams() { }

    public BirthDeathModelParams(double birthRate,
                           double deathRate,
                           double samplingRate,
                           double presentSamplingProb,
                           double originTime) {
        setInputValue(birthRateInput, birthRate);
        setInputValue(deathRateInput, deathRate);
        setInputValue(samplingRateInput, samplingRate);
        setInputValue(presentSamplingProbInput, presentSamplingProb);
        setInputValue(originTimeInput, originTime);
    }

    protected void setInputValue(Input input, double value) {
        input.setValue(new RealParameter(String.valueOf(value)), this);
    }

    @Override
    public void initAndValidate() { }

    public double getBirthRate() {
        return birthRateInput.get().getValue();
    }

    public double getDeathRate() {
        return deathRateInput.get().getValue();
    }

    public double getSamplingRate() {
        return samplingRateInput.get().getValue();
    }

    public double getPresentSamplingProb() {
        return presentSamplingProbInput.get().getValue();
    }

    public double getOriginTime() {
        return originTimeInput.get().getValue();
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
