package com.lifs.tools.ilsceramideringtrial.parser;

import com.lifs.tools.ilsceramideringtrial.model.CeramideConcentrationDataset;
import com.lifs.tools.ilsceramideringtrial.model.CeramideMeasurement;
import com.lifs.tools.ilsceramideringtrial.model.ConcentrationStats;
import com.lifs.tools.ilsceramideringtrial.model.Metadata;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * Parser for ceramide concentration CSV files.
 * Converts CSV data into CeramideConcentrationDataset model objects.
 */
public class CeramideCsvParser {

    // CSV column indices based on the file structure
    private static final int COL_MATRIX = 0;
    private static final int COL_CERAMIDE = 1;
    private static final int COL_N_ALL = 2;
    private static final int COL_MEAN_ALL = 3;
    private static final int COL_SD_ALL = 4;
    private static final int COL_MEDIAN_ALL = 5;
    private static final int COL_CV_ALL = 6;
    private static final int COL_RCV_ALL = 7;
    private static final int COL_N_FILT = 8;
    private static final int COL_MEAN_FILT = 9;
    private static final int COL_SD_FILT = 10;
    private static final int COL_MEDIAN_FILT = 11;
    private static final int COL_CV_FILT = 12;
    private static final int COL_RCV_FILT = 13;

    /**
     * Parses a CSV input stream into a CeramideConcentrationDataset.
     * 
     * @param inputStream The CSV input stream
     * @param metadata The metadata to associate with the dataset
     * @return CeramideConcentrationDataset containing all parsed data
     * @throws IOException If an error occurs during reading
     * @throws CsvParseException If the CSV format is invalid
     */
    public CeramideConcentrationDataset parse(InputStream inputStream, Metadata metadata) 
            throws IOException, CsvParseException {
        
        List<CeramideMeasurement> measurements = new ArrayList<>();
        
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
                
                if (columns.length < 14) {
                    throw new CsvParseException("Invalid CSV line: expected 14 columns, got " + columns.length);
                }
                
                // Extract values
                String matrix = extractValue(columns[COL_MATRIX]);
                String ceramide = extractValue(columns[COL_CERAMIDE]);
                
                // Parse all labs stats
                ConcentrationStats allStats = parseStats(
                    columns[COL_N_ALL],
                    columns[COL_MEAN_ALL],
                    columns[COL_SD_ALL],
                    columns[COL_MEDIAN_ALL],
                    columns[COL_CV_ALL],
                    columns[COL_RCV_ALL]
                );
                
                // Parse filtered stats
                ConcentrationStats filteredStats = parseStats(
                    columns[COL_N_FILT],
                    columns[COL_MEAN_FILT],
                    columns[COL_SD_FILT],
                    columns[COL_MEDIAN_FILT],
                    columns[COL_CV_FILT],
                    columns[COL_RCV_FILT]
                );
                
                // Create measurement
                CeramideMeasurement measurement = new CeramideMeasurement(
                    matrix, ceramide, allStats, filteredStats
                );
                
                measurements.add(measurement);
            }
        }
        
        return new CeramideConcentrationDataset(metadata, measurements);
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
     * Parses statistical values from CSV columns.
     * 
     * @param nCol Number of measurements column
     * @param meanCol Mean concentration column
     * @param sdCol Standard deviation column
     * @param medianCol Median column
     * @param cvCol Coefficient of Variation column
     * @param rcvCol Reference Change Value column
     * @return ConcentrationStats object
     * @throws CsvParseException If any value cannot be parsed
     */
    private ConcentrationStats parseStats(String nCol, String meanCol, String sdCol, 
                                          String medianCol, String cvCol, String rcvCol) 
            throws CsvParseException {
        
        try {
            int n = Integer.parseInt(extractValue(nCol));
            double mean = Double.parseDouble(extractValue(meanCol));
            double sd = Double.parseDouble(extractValue(sdCol));
            double median = Double.parseDouble(extractValue(medianCol));
            double cv = Double.parseDouble(extractValue(cvCol));
            double rcv = Double.parseDouble(extractValue(rcvCol));
            
            return new ConcentrationStats(n, mean, sd, median, cv, rcv);
        } catch (NumberFormatException e) {
            throw new CsvParseException("Failed to parse statistical values: " + e.getMessage(), e);
        }
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
