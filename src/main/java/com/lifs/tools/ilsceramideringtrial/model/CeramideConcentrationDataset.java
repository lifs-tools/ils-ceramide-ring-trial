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
import com.fasterxml.jackson.annotation.JsonProperty;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.AccessLevel;
import lombok.Builder;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;

import java.util.List;

/**
 * Complete dataset model for ceramide concentration measurements.
 * This is the root document that contains metadata and all measurement data.
 * Designed for ArangoDB Spring Data persistence.
 *
 * @author nils.hoffmann
 */
@Schema
@Document(value = "ceramideConcentrationDatasets", keyType = KeyType.autoincrement, keyIncrement = 1, replicationFactor = 2)
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = true)
public class CeramideConcentrationDataset extends ArangoBaseEntity {

    @JsonProperty("metadata")
    private Metadata metadata;

    @JsonProperty("data")
    private List<CeramideMeasurement> data;

    // Constructor
    @Builder
    public CeramideConcentrationDataset(Metadata metadata, List<CeramideMeasurement> data,
                                          String transactionUuid, Visibility visibility) {
        super(transactionUuid, visibility);
        this.metadata = metadata;
        this.data = data;
    }
}

