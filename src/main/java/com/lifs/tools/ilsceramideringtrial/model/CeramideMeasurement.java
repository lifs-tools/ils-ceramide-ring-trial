package com.lifs.tools.ilsceramideringtrial.model;

import com.arangodb.springframework.annotation.Field;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import org.springframework.data.annotation.Id;

/**
 * Represents a single ceramide concentration measurement for a specific matrix.
 * Contains both all-labs data and filtered (outliers removed) data.
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
public class CeramideMeasurement {

    @Id
    private String id;

    @Field("matrix")
    private String matrix;

    @Field("ceramide")
    private String ceramide;

    @Field("all")
    private ConcentrationStats all;

    @Field("filtered")
    private ConcentrationStats filtered;

    public CeramideMeasurement(String matrix, String ceramide, ConcentrationStats all, ConcentrationStats filtered) {
        this.matrix = matrix;
        this.ceramide = ceramide;
        this.all = all;
        this.filtered = filtered;
    }
}
