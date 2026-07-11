/*
 * Copyright 2026 The ILS Ceramide Ring Trial Developers.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.lifs.tools.ilsceramideringtrial.model;

import com.arangodb.entity.KeyType;
import com.arangodb.springframework.annotation.Document;
import com.arangodb.springframework.annotation.PersistentIndexed;
import com.fasterxml.jackson.annotation.JsonProperty;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.AccessLevel;
import lombok.Builder;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;

/**
 * Represents a single lipid concentration measurement for a specific matrix.
 * Contains both all-labs data and filtered (outliers removed) data.
 *
 * @author nils.hoffmann
 */
@Schema
@Document(value = "lipidMeasurements", keyType = KeyType.autoincrement, keyIncrement = 1, replicationFactor = 2)
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = true)
public class LipidMeasurement extends ArangoBaseEntity {

    @PersistentIndexed
    @JsonProperty("matrix")
    private String matrix;

    @PersistentIndexed
    @JsonProperty("lipid")
    private String lipid;

    @JsonProperty("all")
    private ConcentrationStats all;

    @JsonProperty("filtered")
    private ConcentrationStats filtered;

    // Constructor
    @Builder
    public LipidMeasurement(String matrix, String lipid, ConcentrationStats all, 
                            ConcentrationStats filtered, String transactionUuid, Visibility visibility) {
        super(transactionUuid, visibility);
        this.matrix = matrix;
        this.lipid = lipid;
        this.all = all;
        this.filtered = filtered;
    }
}

