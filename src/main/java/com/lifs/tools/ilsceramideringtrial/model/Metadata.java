package com.lifs.tools.ilsceramideringtrial.model;

import com.arangodb.springframework.annotation.Field;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import org.springframework.data.annotation.Id;

import java.util.List;

/**
 * Metadata for the ceramide concentration dataset.
 * Represents publication and dataset information.
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
public class Metadata {

    @Id
    private String id;

    @Field("title")
    private String title;

    @Field("authors")
    private List<String> authors;

    @Field("year")
    private int year;

    @Field("doi")
    private String doi;

    @Field("description")
    private String description;
}
