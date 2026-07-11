package com.lifs.tools.ilsceramideringtrial.model;

import com.arangodb.springframework.annotation.Document;
import com.arangodb.springframework.annotation.Field;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import org.springframework.data.annotation.Id;

import java.util.List;

/**
 * Complete dataset model for ceramide concentration measurements.
 * This is the root document that contains metadata and all measurement data.
 * Designed for ArangoDB Spring Data persistence.
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
@Document(collection = "ceramide_concentration_datasets")
public class CeramideConcentrationDataset {

    @Id
    private String id;

    @Field("metadata")
    private Metadata metadata;

    @Field("data")
    private List<CeramideMeasurement> data;

    public CeramideConcentrationDataset(Metadata metadata, List<CeramideMeasurement> data) {
        this.metadata = metadata;
        this.data = data;
    }
}
