package com.lifs.tools.ilsceramideringtrial.model;

import com.arangodb.springframework.annotation.Field;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;

/**
 * Statistical values for ceramide concentrations.
 * Contains n (count), mean, standard deviation, median, CV, and RCV.
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
public class ConcentrationStats {

    @Field("n")
    private int n;

    @Field("mean")
    private double mean;

    @Field("sd")
    private double sd;

    @Field("median")
    private double median;

    @Field("cv")
    private double cv;

    @Field("rcv")
    private double rcv;

    public ConcentrationStats(int n, double mean, double sd, double median, double cv, double rcv) {
        this.n = n;
        this.mean = mean;
        this.sd = sd;
        this.median = median;
        this.cv = cv;
        this.rcv = rcv;
    }
}
