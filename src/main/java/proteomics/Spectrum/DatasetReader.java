package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import ProteomicsLibrary.MassTool;
import static ProteomicsLibrary.Utilities.*;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.Statement;
import java.util.*;

public class DatasetReader {

    private static final Logger logger = LoggerFactory.getLogger(DatasetReader.class);
    public static final int topN = 15;

    private int usefulSpectraNum = 0;

    public DatasetReader(JMzReader[] spectraParserArray, double ms1Tolerance, int ms1ToleranceUnit, MassTool massTool, String ext, Set<Integer> msLevelSet, String sqlPath, Map<Integer, String> fileIdNameMap) throws Exception {

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        sqlStatement.executeUpdate("PRAGMA journal_mode=WAL");
        sqlStatement.executeUpdate("DROP TABLE IF EXISTS spectraTable");
        sqlStatement.executeUpdate("CREATE TABLE spectraTable (scanNum INTEGER NOT NULL, scanName TEXT PRIMARY KEY, precursorCharge INTEGER NOT NULL, precursorMass REAL NOT NULL, mgfTitle TEXT NOT NULL)");
        sqlStatement.close();

        PreparedStatement sqlPrepareStatement = sqlConnection.prepareStatement("INSERT INTO spectraTable (scanNum, scanName, precursorCharge, precursorMass, mgfTitle) VALUES (?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);

        for (int i = 0; i < spectraParserArray.length; i++) {
            JMzReader spectraParser = spectraParserArray[i];
            Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();

            while (spectrumIterator.hasNext()) {
                try {
                    Spectrum spectrum = spectrumIterator.next();
                    if (ext.toLowerCase().contentEquals("mzxml")) {
                        if (!msLevelSet.contains(spectrum.getMsLevel())) {
                            continue;
                        }
                    }

                    int scanNum;
                    double precursorMz = spectrum.getPrecursorMZ();
                    int precursorCharge = -1;
                    double precursorMass;
                    String mgfTitle = "";
                    mgfTitle = ((Ms2Query) spectrum).getTitle();
                    scanNum = getScanNum(mgfTitle);

                    if (spectrum.getPeakList().size() < 5) {
                        continue;
                    }

                    if (spectrum.getPrecursorCharge() == null) {
                        continue;
                    } else {
                        precursorCharge = spectrum.getPrecursorCharge();
                        precursorMass = precursorMz * precursorCharge - precursorCharge * MassTool.PROTON;
                    }
                    sqlPrepareStatement.setInt(1, scanNum);
                    sqlPrepareStatement.setString(2, fileIdNameMap.get(i)+"."+spectrum.getId()+"."+scanNum); //fileName.scanId.scanNum  scanName
                    sqlPrepareStatement.setInt(3, precursorCharge);
                    sqlPrepareStatement.setDouble(4, precursorMass);
                    sqlPrepareStatement.setString(5, mgfTitle);
                    sqlPrepareStatement.executeUpdate();
                    ++usefulSpectraNum;
                } catch (RuntimeException ex) {
                    logger.error(ex.toString());
                }
            }
        }

        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        sqlPrepareStatement.close();
        sqlConnection.close();
        logger.info("Useful MS/MS spectra number: {}.", usefulSpectraNum);
    }
    public int getUsefulSpectraNum() {
        return usefulSpectraNum;
    }
}
