#pragma once

// Shasta
#include <invalid.hpp>

// Standard library.
#include <cstdint.hpp>
#include <iosfwd.hpp>
#include <vector.hpp>



// Likelihood ratio test of the tangle matrix (G test).
// https://en.wikipedia.org/wiki/G-test
namespace shasta {
    class GTest;
}



class shasta::GTest {
public:
    GTest(const vector< vector<uint64_t> >& tangleMatrix, double epsilon);
    GTest(const vector< vector<double> >& tangleMatrix, double epsilon);
    bool success = false;

    void writeHtml(ostream&) const;

    class Hypothesis {
    public:
        vector< vector<bool> > connectivityMatrix;
        double G = invalid<double>;

        Hypothesis(
            const vector< vector<bool> >& connectivityMatrix,
            double G) :
            connectivityMatrix(connectivityMatrix),
            G(G)
            {}

        Hypothesis() {}

        // Return true if there is a single exit for each entrance.
        bool isForwardInjective() const;

        // Return true if there is a single entrance for exit entrance.
        bool isBackwardInjective() const;

        // Sort by G.
        bool operator<(const Hypothesis& that) const {
            return G < that.G;
        }
    };
    vector<Hypothesis> hypotheses;

private:
    void run(const vector< vector<double> >& tangleMatrix, double epsilon);

};
