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

import com.arangodb.springframework.annotation.ArangoId;
import com.arangodb.springframework.annotation.Rev;
import com.fasterxml.jackson.annotation.JsonProperty;
import lombok.AccessLevel;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;
import lombok.experimental.SuperBuilder;
import org.springframework.data.annotation.CreatedBy;
import org.springframework.data.annotation.CreatedDate;
import org.springframework.data.annotation.Id;
import org.springframework.data.annotation.LastModifiedBy;
import org.springframework.data.annotation.LastModifiedDate;

import javax.validation.constraints.NotNull;
import java.util.Date;

/**
 * Base entity class for ArangoDB documents.
 * Provides common fields for all persisted entities.
 *
 * @author nils.hoffmann
 */
@Data
@AllArgsConstructor
@NoArgsConstructor(access = AccessLevel.PROTECTED)
@EqualsAndHashCode(onlyExplicitlyIncluded = true)
@SuperBuilder
public abstract class ArangoBaseEntity {

    @Id
    private String id;

    @ArangoId
    private String arangoId;

    @JsonProperty("revision")
    @Rev
    private String revision;

    @NotNull(message = "transactionUuid must not be null!")
    private String transactionUuid;

    @JsonProperty("dateCreated")
    @CreatedDate
    private Date dateCreated;

    @JsonProperty("dateLastModified")
    @LastModifiedDate
    private Date dateLastModified;

    @JsonProperty("createdBy")
    @CreatedBy
    private String createdBy;

    @JsonProperty("updatedBy")
    @LastModifiedBy
    private String updatedBy;

    @Builder.Default
    private Visibility visibility = Visibility.PRIVATE;

    public ArangoBaseEntity(String transactionUuid, Visibility visibility) {
        if (transactionUuid == null) {
            throw new IllegalArgumentException("Argument transactionUuid must not be null!");
        }
        if (visibility == null) {
            throw new IllegalArgumentException("Argument visibility must not be null!");
        }
        this.transactionUuid = transactionUuid;
        this.visibility = visibility;
    }
}
