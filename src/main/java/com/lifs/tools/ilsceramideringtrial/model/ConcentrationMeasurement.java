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
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;

/**
 * Represents a single concentration measurement from a lab in a ring trial.
 * Contains the measured quantity and metadata about the measurement group.
 *
 * @author nils.hoffmann
 */
@Schema
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@AllArgsConstructor
@Builder
public class ConcentrationMeasurement {

    @JsonProperty("sample_matrix")
    private String sampleMatrix;

    @JsonProperty("lipid")
    private String lipid;

    @JsonProperty("source_id")
    private String sourceId;

    @JsonProperty("quantity")
    private double quantity;

    @JsonProperty("group")
    private String group;
}
