package hobd.simulator;

public class HOBDParams {
    public double birthRate, burstRate, burstSize, deathRate;
    public double psiSampleRate, rhoSampleProb;
    public double T;

    public HOBDParams(double birthRate,
                                   double deathRate,
                                   double burstRate,
                                   double burstSize,
                                   double psiSampleRate,
                                   double rhoSampleProb,
                                   double T) {

        this.birthRate = birthRate;
        this.burstRate = burstRate;
        this.burstSize = burstSize;
        this.deathRate = deathRate;
        this.psiSampleRate = psiSampleRate;
        this.rhoSampleProb = rhoSampleProb;
        this.T = T;
    }
}
