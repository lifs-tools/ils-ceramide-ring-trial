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
 * Represents a single ceramide concentration measurement for a specific matrix.
 * Contains both all-labs data and filtered (outliers removed) data.
 *
 * @author nils.hoffmann
 */
@Schema
@Document(value = "ceramideMeasurements", keyType = KeyType.autoincrement, keyIncrement = 1, replicationFactor = 2)
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = true)
public class CeramideMeasurement extends ArangoBaseEntity {

    @PersistentIndexed
    @JsonProperty("matrix")
    private String matrix;

    @PersistentIndexed
    @JsonProperty("ceramide")
    private String ceramide;

    @JsonProperty("all")
    private ConcentrationStats all;

    @JsonProperty("filtered")
    private ConcentrationStats filtered;

    // Constructor
    @Builder
    public CeramideMeasurement(String matrix, String ceramide, ConcentrationStats all, 
                                ConcentrationStats filtered, String transactionUuid, Visibility visibility) {
        super(transactionUuid, visibility);
        this.matrix = matrix;
        this.ceramide = ceramide;
        this.all = all;
        this.filtered = filtered;
    }
}

