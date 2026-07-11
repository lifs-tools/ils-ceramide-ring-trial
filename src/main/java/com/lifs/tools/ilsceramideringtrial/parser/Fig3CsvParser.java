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
package com.lifs.tools.ilsceramideringtrial.parser;

import com.lifs.tools.ilsceramideringtrial.model.ConcentrationMeasurement;
import com.lifs.tools.ilsceramideringtrial.model.ConcentrationStats;
import com.lifs.tools.ilsceramideringtrial.model.PublicationMetadata;
import com.lifs.tools.ilsceramideringtrial.model.RingTrialMeasurement;
import com.lifs.tools.ilsceramideringtrial.model.RingTrialResult;
import com.lifs.tools.ilsceramideringtrial.model.Visibility;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

/**
 * Parser for Fig3 source data CSV files.
 * Converts CSV data into RingTrialResult model objects with individual
 * lab measurements grouped by sample matrix and lipid.
 *
 * @author nils.hoffmann
 */
public class Fig3CsvParser {

    // CSV column indices based on the Fig3 file structure
    private static final int COL_SAMPLE_TYPE = 0;
    private static final int COL_CERAMIDE_NAME = 1;
    private static final int COL_LAB_ID = 2;
    private static final int COL_LAB_ID_GROUP = 3;
    private static final int COL_LAB_NUM = 4;
    private static final int COL_METHOD_NO = 5;
    private static final int COL_PROTOCOL = 6;
    private static final int COL_LC = 7;
    private static final int COL_MASS_ANALYZER_RESOLUTION = 8;
    private static final int COL_MASS_ANALYZER_TYPE = 9;
    private static final int COL_C_ADJ_MEAN = 14;

    /**
     * Parses a CSV input stream into a RingTrialResult.
     * 
     * @param inputStream The CSV input stream
     * @param publicationMetadata The publication metadata to associate with the dataset
     * @return RingTrialResult containing all parsed data
     * @throws IOException If an error occurs during reading
     * @throws CsvParseException If the CSV format is invalid
     */
    public RingTrialResult parse(InputStream inputStream, PublicationMetadata publicationMetadata) 
            throws IOException, CsvParseException {
        
        // Map to group measurements by (sampleMatrix, lipid) -> RingTrialMeasurement
        Map<String, Map<String, RingTrialMeasurement.Builder>> measurementMap = new HashMap<>();
        
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(inputStream, StandardCharsets.UTF_8))) {
            
            String line;
            boolean isFirstLine = true;
            
            while ((line = reader.readLine()) != null) {
                // Skip empty lines
                if (line.trim().isEmpty()) {
                    continue;
                }
                
                // Skip header line (first line)
                if (isFirstLine) {
                    isFirstLine = false;
                    continue;
                }
                
                // Parse the CSV line
                String[] columns = parseCsvLine(line);
                
                if (columns.length < 15) {
                    throw new CsvParseException("Invalid CSV line: expected at least 15 columns, got " + columns.length);
                }
                
                // Extract values
                String sampleMatrix = extractValue(columns[COL_SAMPLE_TYPE]);
                String lipid = extractValue(columns[COL_CERAMIDE_NAME]);
                String sourceId = extractValue(columns[COL_LAB_NUM]);
                double quantity = parseDouble(columns[COL_C_ADJ_MEAN]);
                
                // Build group string from Protocol, LC, MassAnalyzerResolution, MassAnalyzerType
                String group = buildGroupString(
                    columns[COL_PROTOCOL],
                    columns[COL_LC],
                    columns[COL_MASS_ANALYZER_RESOLUTION],
                    columns[COL_MASS_ANALYZER_TYPE]
                );
                
                // Create concentration measurement
                ConcentrationMeasurement concentrationMeasurement = ConcentrationMeasurement.builder()
                    .sampleMatrix(sampleMatrix)
                    .lipid(lipid)
                    .sourceId(sourceId)
                    .quantity(quantity)
                    .group(group)
                    .build();
                
                // Get or create the RingTrialMeasurement builder for this (sampleMatrix, lipid) combination
                RingTrialMeasurement.Builder rtMeasurementBuilder = measurementMap
                    .computeIfAbsent(sampleMatrix, k -> new HashMap<>())
                    .computeIfAbsent(lipid, k -> RingTrialMeasurement.builder()
                        .sampleMatrix(sampleMatrix)
                        .lipid(lipid)
                        .all(new ConcentrationStats(0, 0, 0, 0, 0, 0))  // Placeholder, will be calculated
                        .filtered(new ConcentrationStats(0, 0, 0, 0, 0, 0))  // Placeholder, will be calculated
                        .measurements(new ArrayList<>())
                        .transactionUuid(UUID.randomUUID().toString())
                        .visibility(Visibility.PUBLIC));
                
                // Add the concentration measurement to the list
                rtMeasurementBuilder.measurements.add(concentrationMeasurement);
            }
        }
        
        // Build the final list of RingTrialMeasurements
        List<RingTrialMeasurement> measurements = new ArrayList<>();
        for (Map<String, RingTrialMeasurement.Builder> lipidMap : measurementMap.values()) {
            for (RingTrialMeasurement.Builder builder : lipidMap.values()) {
                // Calculate statistics from the measurements
                ConcentrationStats allStats = calculateStats(builder.measurements);
                builder.all(allStats);
                builder.filtered(allStats);  // For now, use same stats for filtered
                
                measurements.add(builder.build());
            }
        }
        
        // Build the final result
        return RingTrialResult.builder()
            .publicationMetadata(publicationMetadata)
            .data(measurements)
            .transactionUuid(UUID.randomUUID().toString())
            .visibility(Visibility.PUBLIC)
            .build();
    }

    /**
     * Builds a group string from the protocol, LC, mass analyzer resolution, and mass analyzer type.
     *
     * @param protocol The protocol
     * @param lc The LC type
     * @param massAnalyzerResolution The mass analyzer resolution
     * @param massAnalyzerType The mass analyzer type
     * @return Combined group string
     */
    private String buildGroupString(String protocol, String lc, String massAnalyzerResolution, String massAnalyzerType) {
        return String.format("%s | %s | %s | %s",
            extractValue(protocol),
            extractValue(lc),
            extractValue(massAnalyzerResolution),
            extractValue(massAnalyzerType));
    }

    /**
     * Calculates statistics from a list of concentration measurements.
     *
     * @param measurements The list of measurements
     * @return ConcentrationStats with calculated values
     */
    private ConcentrationStats calculateStats(List<ConcentrationMeasurement> measurements) {
        if (measurements == null || measurements.isEmpty()) {
            return new ConcentrationStats(0, 0, 0, 0, 0, 0);
        }
        
        int n = measurements.size();
        
        // Calculate mean
        double sum = 0;
        for (ConcentrationMeasurement m : measurements) {
            sum += m.getQuantity();
        }
        double mean = sum / n;
        
        // Calculate standard deviation
        double sumSquaredDiff = 0;
        for (ConcentrationMeasurement m : measurements) {
            double diff = m.getQuantity() - mean;
            sumSquaredDiff += diff * diff;
        }
        double sd = Math.sqrt(sumSquaredDiff / n);
        
        // Calculate median
        List<Double> sortedQuantities = new ArrayList<>();
        for (ConcentrationMeasurement m : measurements) {
            sortedQuantities.add(m.getQuantity());
        }
        sortedQuantities.sort(Double::compareTo);
        double median = sortedQuantities.get(sortedQuantities.size() / 2);
        
        // Calculate CV (Coefficient of Variation)
        double cv = (mean != 0) ? (sd / mean) * 100 : 0;
        
        // RCV is not calculated here, set to 0
        double rcv = 0;
        
        return new ConcentrationStats(n, mean, sd, median, cv, rcv);
    }

    /**
     * Parses a double value from a CSV column.
     *
     * @param value The column value
     * @return Parsed double value, or 0 if NA or cannot parse
     */
    private double parseDouble(String value) {
        String cleaned = extractValue(value);
        if (cleaned == null || cleaned.equalsIgnoreCase("NA") || cleaned.isEmpty()) {
            return 0.0;
        }
        try {
            return Double.parseDouble(cleaned);
        } catch (NumberFormatException e) {
            return 0.0;
        }
    }

    /**
     * Parses a CSV line handling quoted values.
     * 
     * @param line The CSV line
     * @return Array of column values
     */
    private String[] parseCsvLine(String line) {
        List<String> columns = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;
        
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                columns.add(current.toString());
                current = new StringBuilder();
            } else {
                current.append(c);
            }
        }
        
        // Add the last column
        columns.add(current.toString());
        
        return columns.toArray(new String[0]);
    }

    /**
     * Extracts and cleans a value from CSV column.
     * Removes surrounding quotes and whitespace.
     * 
     * @param value The raw CSV column value
     * @return Cleaned value
     */
    private String extractValue(String value) {
        if (value == null) {
            return null;
        }
        
        String cleaned = value.trim();
        
        // Remove surrounding quotes if present
        if (cleaned.startsWith("\"") && cleaned.endsWith("\"")) {
            cleaned = cleaned.substring(1, cleaned.length() - 1);
        }
        
        return cleaned;
    }

    /**
     * Exception thrown when CSV parsing fails.
     */
    public static class CsvParseException extends Exception {
        public CsvParseException(String message) {
            super(message);
        }

        public CsvParseException(String message, Throwable cause) {
            super(message, cause);
        }
    }
}
