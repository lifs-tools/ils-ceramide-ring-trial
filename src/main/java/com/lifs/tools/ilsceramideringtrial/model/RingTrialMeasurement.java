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

import com.fasterxml.jackson.annotation.JsonProperty;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.AccessLevel;
import lombok.Builder;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;

import java.util.List;

/**
 * Represents a lipid measurement in a ring trial with individual lab measurements.
 * Extends LipidMeasurement to add the list of concentration measurements from
 * individual labs.
 *
 * @author nils.hoffmann
 */
@Schema
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = true)
public class RingTrialMeasurement extends LipidMeasurement {

    @JsonProperty("measurements")
    private List<ConcentrationMeasurement> measurements;

    // Constructor
    @Builder
    public RingTrialMeasurement(String sampleMatrix, String lipid, ConcentrationStats all, 
                                ConcentrationStats filtered, List<ConcentrationMeasurement> measurements,
                                String transactionUuid, Visibility visibility) {
        super(sampleMatrix, lipid, all, filtered, transactionUuid, visibility);
        this.measurements = measurements;
    }
}
