package hobd;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

public class HOBDModelParams extends CalculationNode {

    public Input<RealParameter> birthRateInput = new Input<>("birthRate",
            "Birth rate", Input.Validate.REQUIRED);

    public Input<RealParameter> deathRateInput = new Input<>("deathRate",
            "Death rate", Input.Validate.REQUIRED);

    public Input<RealParameter> burstRateInput = new Input<>("burstRate",
            "Burst rate", Input.Validate.REQUIRED);

    public Input<RealParameter> burstSizeInput = new Input<>("burstSize",
            "Mean burst rate", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingRateInput = new Input<>("samplingRate",
            "Sampling rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> presentSamplingProbInput = new Input<>("presentSamplingProb",
            "Present-day sampling prob.", Input.Validate.REQUIRED);

    public Input<RealParameter> originTimeInput = new Input<>("originTime",
            "Time between present and origin of process.",
            Input.Validate.REQUIRED);

    public HOBDModelParams() { }

    public HOBDModelParams(double birthRate,
                           double deathRate,
                           double burstRate,
                           double burstSize,
                           double samplingRate,
                           double presentSamplingProb,
                           double originTime) {
        setInputValue(birthRateInput, birthRate);
        setInputValue(deathRateInput, deathRate);
        setInputValue(burstRateInput, burstRate);
        setInputValue(burstSizeInput, burstSize);
        setInputValue(samplingRateInput, samplingRate);
        setInputValue(presentSamplingProbInput, presentSamplingProb);
        setInputValue(originTimeInput, originTime);
    }

    private void setInputValue(Input input, double value) {
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

    public double getBurstRate() {
        return burstRateInput.get().getValue();
    }

    public double getBurstSize() {
        return burstSizeInput.get().getValue();
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
