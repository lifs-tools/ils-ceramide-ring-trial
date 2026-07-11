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
 * Represents ring trial results for lipid concentration measurements.
 * Similar to LipidConcentrationDataset but with data being a list of RingTrialMeasurement
 * objects that include individual lab measurements.
 *
 * @author nils.hoffmann
 */
@Schema
@Document(value = "ringTrialResults", keyType = KeyType.autoincrement, keyIncrement = 1, replicationFactor = 2)
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = false)
public class RingTrialResult extends ArangoBaseEntity {

    @JsonProperty("publicationMetadata")
    private PublicationMetadata publicationMetadata;

    @JsonProperty("data")
    private List<RingTrialMeasurement> data;

    // Constructor
    @Builder
    public RingTrialResult(PublicationMetadata publicationMetadata, 
                          List<RingTrialMeasurement> data,
                          String transactionUuid, Visibility visibility) {
        super(transactionUuid, visibility);
        this.publicationMetadata = publicationMetadata;
        this.data = data;
    }
}
