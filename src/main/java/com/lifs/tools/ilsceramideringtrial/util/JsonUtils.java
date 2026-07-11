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
package com.lifs.tools.ilsceramideringtrial.util;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.lifs.tools.ilsceramideringtrial.model.CeramideConcentrationDataset;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Utility class for JSON serialization and deserialization.
 * Uses Jackson for JSON processing.
 */
public class JsonUtils {

    private static final ObjectMapper objectMapper = createObjectMapper();

    private static ObjectMapper createObjectMapper() {
        ObjectMapper mapper = new ObjectMapper();
        mapper.enable(SerializationFeature.INDENT_OUTPUT);
        mapper.findAndRegisterModules();
        return mapper;
    }

    /**
     * Serializes an object to JSON string.
     * 
     * @param object The object to serialize
     * @return JSON string representation
     * @throws JsonProcessingException If serialization fails
     */
    public static String toJson(Object object) throws JsonProcessingException {
        return objectMapper.writeValueAsString(object);
    }

    /**
     * Deserializes JSON string to an object.
     * 
     * @param json JSON string
     * @param valueType The target class type
     * @param <T> The target type
     * @return Deserialized object
     * @throws JsonProcessingException If deserialization fails
     */
    public static <T> T fromJson(String json, Class<T> valueType) throws JsonProcessingException {
        return objectMapper.readValue(json, valueType);
    }

    /**
     * Serializes an object to a file.
     * 
     * @param object The object to serialize
     * @param file The target file
     * @throws IOException If writing fails
     */
    public static void toFile(Object object, File file) throws IOException {
        objectMapper.writeValue(file, object);
    }

    /**
     * Deserializes a file to an object.
     * 
     * @param file The JSON file
     * @param valueType The target class type
     * @param <T> The target type
     * @return Deserialized object
     * @throws IOException If reading fails
     */
    public static <T> T fromFile(File file, Class<T> valueType) throws IOException {
        return objectMapper.readValue(file, valueType);
    }

    /**
     * Serializes an object to an output stream.
     * 
     * @param object The object to serialize
     * @param outputStream The target output stream
     * @throws IOException If writing fails
     */
    public static void toStream(Object object, OutputStream outputStream) throws IOException {
        objectMapper.writeValue(outputStream, object);
    }

    /**
     * Deserializes an input stream to an object.
     * 
     * @param inputStream The JSON input stream
     * @param valueType The target class type
     * @param <T> The target type
     * @return Deserialized object
     * @throws IOException If reading fails
     */
    public static <T> T fromStream(InputStream inputStream, Class<T> valueType) throws IOException {
        return objectMapper.readValue(inputStream, valueType);
    }

    /**
     * Validates JSON against the schema.
     * Note: This is a placeholder for actual schema validation.
     * For full schema validation, use a library like NetworkNT/schema-validator.
     * 
     * @param json JSON string to validate
     * @return true if valid, false otherwise
     */
    public static boolean validateAgainstSchema(String json) {
        try {
            // Basic validation - try to parse as CeramideConcentrationDataset
            objectMapper.readValue(json, CeramideConcentrationDataset.class);
            return true;
        } catch (JsonProcessingException e) {
            return false;
        }
    }

    /**
     * Pretty prints JSON string.
     * 
     * @param json JSON string to format
     * @return Formatted JSON string
     * @throws JsonProcessingException If parsing fails
     */
    public static String prettyPrint(String json) throws JsonProcessingException {
        Object value = objectMapper.readValue(json, Object.class);
        return objectMapper.writerWithDefaultPrettyPrinter().writeValueAsString(value);
    }
}
