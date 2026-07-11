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

import com.arangodb.springframework.annotation.PersistentIndexed;
import com.fasterxml.jackson.annotation.JsonProperty;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.AccessLevel;
import lombok.Builder;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;

import java.util.List;

/**
 * Metadata for the ceramide concentration dataset.
 * Represents publication and dataset information.
 *
 * @author nils.hoffmann
 */
@Schema
@Data
@NoArgsConstructor(access = AccessLevel.PRIVATE)
@EqualsAndHashCode(onlyExplicitlyIncluded = true, callSuper = true)
public class Metadata extends ArangoBaseEntity {

    @PersistentIndexed
    @JsonProperty("title")
    private String title;

    @JsonProperty("authors")
    private List<String> authors;

    @PersistentIndexed
    @JsonProperty("year")
    private int year;

    @PersistentIndexed
    @JsonProperty("doi")
    private String doi;

    @JsonProperty("description")
    private String description;

    // Constructor
    @Builder
    public Metadata(String title, List<String> authors, int year, String doi, String description,
                    String transactionUuid, Visibility visibility) {
        super(transactionUuid, visibility);
        this.title = title;
        this.authors = authors;
        this.year = year;
        this.doi = doi;
        this.description = description;
    }
}
