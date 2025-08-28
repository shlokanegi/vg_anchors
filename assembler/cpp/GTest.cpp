// Shasta.
#include "GTest.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

// Standard library.
#include <algorithm.hpp>
#include <cmath>
#include <iomanip>
#include <iostream.hpp>
#include <numeric>



GTest::GTest(const vector< vector<double> >& tangleMatrix, double epsilon)
{
    run(tangleMatrix, epsilon);
}



GTest::GTest(
    const vector< vector<uint64_t> >& tangleMatrixInteger,
    double epsilon)
{
    // Make a copy of the tangle matrix that stores double instead of uint64_t.
    vector< vector<double> > tangleMatrixDouble;
    for(const vector<uint64_t>& rowInteger: tangleMatrixInteger) {
        vector<double> rowDouble;
        for(const uint64_t valueInteger: rowInteger) {
            rowDouble.push_back(double(valueInteger));
        }
        tangleMatrixDouble.push_back(rowDouble);
    }

    run(tangleMatrixDouble, epsilon);
}



void GTest::run(const vector< vector<double> >& tangleMatrix, double epsilon)
{

    const bool debug = false;

    // Get the number of entrances.
    // This is the number of rows in the tangle matrix.
    const uint64_t entranceCount = tangleMatrix.size();
    SHASTA_ASSERT(entranceCount > 0);

    // Get the number of exits.
    // This is the number of columns in the tangle matrix.
    const uint64_t exitCount = tangleMatrix.front().size();
    SHASTA_ASSERT(exitCount > 0);
    for(uint64_t i=0; i<entranceCount;i++) {
        SHASTA_ASSERT(tangleMatrix[i].size() == exitCount);
    }

    // Limit to a maximum of 16 tangle matrix entries.
    const uint64_t totalTangleMatrixEntryCount = entranceCount * exitCount;
    if(totalTangleMatrixEntryCount > 16) {
        success = false;
        return;
    }

    // Compute common coverage for each entrance and for each exit.
    // Common coverage for an entrance is the sum the tangle matrix row for that entrance.
    // Common coverage for an exit is the sum the tangle matrix column for that exit.
    double totalCommonCoverage = 0;
    vector<double> entranceCommonCoverage(entranceCount, 0UL);
    vector<double> exitCommonCoverage(exitCount, 0UL);
    for(uint64_t i=0; i<entranceCount; i++) {
        for(uint64_t j=0; j<exitCount; j++) {
            const double coverage = tangleMatrix[i][j];
            totalCommonCoverage += coverage;
            entranceCommonCoverage[i] += coverage;
            exitCommonCoverage[j] += coverage;
        }
    }

    // Compute what the tangle matrix would be under entirely random assumptions.
    vector< vector<double> > randomTangleMatrix(entranceCount, vector<double>(exitCount, 0.));
    for(uint64_t i=0; i<entranceCount; i++) {
        for(uint64_t j=0; j<exitCount; j++) {
            randomTangleMatrix[i][j] =
                entranceCommonCoverage[i] *
                exitCommonCoverage[j] /
                totalCommonCoverage;
        }
    }

    // Some matrices used in the loop below.
    vector< vector<double> > idealTangleMatrixA(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrixB(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<double> > expectedTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<bool> > connectivityMatrix(entranceCount, vector<bool>(exitCount));



    // The connectivity matrix contains true for entrance/exit pairs to be connected.
    // Try all N possible connectivity matrices. Each of them can generate a Hypothesis.
    const uint64_t N = 1ULL << totalTangleMatrixEntryCount;
    hypotheses.clear();
    for(uint64_t connectivityInteger=0; connectivityInteger<N; connectivityInteger++) {

        // Use the bits of connectivityInteger to construct the connectivity matrix.
        uint64_t mask = 1;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                connectivityMatrix[iEntrance][iExit] = ((connectivityInteger & mask) != 0);
                mask = mask << 1;
            }
        }

        if(debug) {
            cout << "Trying connectivity matrix:" << endl;
            for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
                for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                    cout << connectivityMatrix[iEntrance][iExit] << " ";
                }
                cout << endl;
            }
        }

        // Compute the tangle matrix we would see assuming this connectivity matrix
        // and no errors. This can be done in two ways:
        // - Equally distributing common coverage at each entrance among all the
        //   exits for which the connectivity matrix is true for that entrance (idealTangleMatrixA).
        // - Equally distributing common coverage at each exit among all the
        //   entrances for which the connectivity matrix is true for that entrance (idealTangleMatrixB).
        // We average the two ways to compute the idealTangleMatrix.

        // Compute idealTangleMatrixA.
        for(uint64_t i=0; i<entranceCount; i++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t j=0; j<exitCount; j++) {
                if(connectivityMatrix[i][j]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t j=0; j<exitCount; j++) {
                    idealTangleMatrixA[i][j] = 0.;
                }
            } else {
                const double value = entranceCommonCoverage[i] / double(nonZeroCount);
                for(uint64_t j=0; j<exitCount; j++) {
                    if(connectivityMatrix[i][j]) {
                        idealTangleMatrixA[i][j] = value;
                    } else {
                        idealTangleMatrixA[i][j] = 0.;
                    }
                }
            }
        }

        // Compute idealTangleMatrixB.
        for(uint64_t j=0; j<exitCount; j++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t i=0; i<entranceCount; i++) {
                if(connectivityMatrix[i][j]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t i=0; i<entranceCount; i++) {
                    idealTangleMatrixB[i][j] = 0.;
                }
            } else {
                const double value = exitCommonCoverage[j] / double(nonZeroCount);
                for(uint64_t i=0; i<entranceCount; i++) {
                    if(connectivityMatrix[i][j]) {
                        idealTangleMatrixB[i][j] = value;
                    } else {
                        idealTangleMatrixB[i][j] = 0.;
                    }
                }
            }
        }

        // Compute the ideal tangle matrix by averaging idealTangleMatrixA and idealTangleMatrixB.
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                idealTangleMatrix[iEntrance][iExit] =  0.5 *
                    (idealTangleMatrixA[iEntrance][iExit] + idealTangleMatrixB[iEntrance][iExit]);
            }
        }

        // Compute the expected tangle matrix given this connectivity matrix.
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                expectedTangleMatrix[i][j] =
                    (1. - epsilon) * idealTangleMatrix [i][j] +
                    epsilon        * randomTangleMatrix[i][j];
            }
        }

        if(debug) {
            cout << "Expected non-ideal tangle matrix for this connectivity matrix:" << endl;
            for(uint64_t i=0; i<entranceCount; i++) {
                for(uint64_t j=0; j<exitCount; j++) {
                    cout << expectedTangleMatrix[i][j] << " ";
                }
                cout << endl;
            }
        }

        // Do a G-test of the observed tangle matrix against the expected one.
        // https://en.wikipedia.org/wiki/G-test
        double G = 0.;
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                const double actualCoverage = tangleMatrix[i][j];
                if(actualCoverage > 0.) {
                    const double expectedCoverage = expectedTangleMatrix[i][j];
                    G += actualCoverage * log10(actualCoverage / expectedCoverage);
                }
            }
        }
        G *= 20.;  // Factor of 2 and convert to decibels.

        if(debug) {
            cout << "G " << G << endl;
        }

        // Store this hypothesis.
        hypotheses.emplace_back(Hypothesis(connectivityMatrix, G));

    }

    sort(hypotheses.begin(), hypotheses.end());
    success = not hypotheses.empty();
}



void GTest::writeHtml(ostream& html) const
{
    html <<
        "<h3>G test</h3>"
        "<table><tr>"
        "<th>Connectivity<br>matrix"
        "<th>G<br>(dB)";

    for(const auto& hypothesis: hypotheses) {
        const vector< vector<bool> >& connectivityMatrix = hypothesis.connectivityMatrix;

        html << "<tr><td style='display: flex; align-items: center; justify-content: center;'>";

        html << "<table>";
        for(uint64_t i=0; i<connectivityMatrix.size(); i++) {
            html << "<tr>";
            for(uint64_t j=0; j<connectivityMatrix[i].size(); j++) {
                html << "<td class=centered>" << hypothesis.connectivityMatrix[i][j];
            }
        }
        html << "</table>";

        html <<
            "<td class=centered>" << std::fixed << std::setprecision(1) << hypothesis.G;

    }

    html << "</table>";
}



// Return true if there is a single exit for each entrance.
bool GTest::Hypothesis::isForwardInjective() const
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        uint64_t count = 0;
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;
}



// Return true if there is a single entrance for exit entrance.
bool GTest::Hypothesis::isBackwardInjective() const
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        uint64_t count = 0;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;

}
