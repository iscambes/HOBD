package hobd;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class HOBDModelParams extends BirthDeathModelParams {

    public Input<RealParameter> burstRateInput = new Input<>("burstRate",
            "Per-individual burst rate", Input.Validate.REQUIRED);

    public Input<RealParameter> meanBurstSizeInput = new Input<>("meanBurstSize",
            "Mean burst size (nr. of new infections - 1)", Input.Validate.REQUIRED);

    public HOBDModelParams() { }

    public HOBDModelParams(double birthRate,
                           double deathRate,
                           double burstRate,
                           double meanBurstSize,
                           double samplingRate,
                           double presentSamplingProb,
                           double originTime) {
        setInputValue(birthRateInput, birthRate);
        setInputValue(deathRateInput, deathRate);
        setInputValue(burstRateInput, burstRate);
        setInputValue(meanBurstSizeInput, meanBurstSize);
        setInputValue(samplingRateInput, samplingRate);
        setInputValue(presentSamplingProbInput, presentSamplingProb);
        setInputValue(originTimeInput, originTime);
    }

    public double getBurstRate() {
        return burstRateInput.get().getValue();
    }

    public double getMeanBurstSize() {
        return meanBurstSizeInput.get().getValue();
    }
}
