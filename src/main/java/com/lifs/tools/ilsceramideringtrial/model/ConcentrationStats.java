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
 * Statistical values for ceramide concentrations.
 * Contains n (count), mean, standard deviation, median, CV, and RCV.
 *
 * @author nils.hoffmann
 */
@Schema
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@AllArgsConstructor
@Builder
public class ConcentrationStats {

    @JsonProperty("n")
    private int n;

    @JsonProperty("mean")
    private double mean;

    @JsonProperty("sd")
    private double sd;

    @JsonProperty("median")
    private double median;

    @JsonProperty("cv")
    private double cv;

    @JsonProperty("rcv")
    private double rcv;
}

