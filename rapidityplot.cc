// PYTHIA script to generate histograms of rapidities of primary hadrons
// produced by fragmentation. Simulation is done for a single q-qbar string
// and only the hadronisation process is considered, with parton shower and
// other effects disabled. The invariant mass of the string can be varied.

// Author: Jade Abidi
// Created: 30/07/2025

#include "Pythia8/Pythia.h"
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

using namespace Pythia8;
using namespace std;

// Function to convert string masses to string with 2 decimal places.
string to_string_2dp(double val) {
  ostringstream out;
  out << fixed << setprecision(2) << val;
  return out.str();
}

int main() {
  // Specify invariant string CMEs to be simulated (in GeV), and number of
  // events to simulate per string CME.
  vector<double> masses = {5, 20, 100};
  vector<string> plotColours = {"steelblue", "seagreen", "indianred"};
  int nEvent = 1000000;

  // Specify id of quark.
  // 1 - down. 2 - up. 3 - strange. 4 - charm. 5 - bottom. 6 - top.
  int qid = 1;

  // Option for massless quarks.
  bool masslessQuarks = true;

  // Initialise matplotlib setup
  HistPlot hpl("rapidityplot");
  hpl.frame("rapidityplot", "Rapidity distributions of primary "
	    "hadrons for differing string energies", "y", "n");

  // Run separately for each invariant string mass.
  for (size_t iMass = 0; iMass < masses.size(); ++iMass) {
    double smass = masses[iMass];

    // Set up generator.
    Pythia pythia;
    Event& event = pythia.event;
    ParticleData& pdt = pythia.particleData;

    // Disable parton shower and hard process since q-qbar will be manually
    // input.
    pythia.readString("ProcessLevel:all = off");

    // Disable hadron decay.
    pythia.readString("HadronLevel:Decay = off");

    // Optional: Disable transverse momentum (enforce 1+1 dimensions)
    pythia.readString("StringPT:sigma = 0");

    // Customise output to be more readable and less cluttered.
    pythia.readString("Next:numberCount = 100000");

    // Initialise.
    cout << "Initialising PYTHIA for q-qbar hadronisation, string mass = " \
	 << smass << endl;
    if (!pythia.init()) return 1;
    
    // Set up histogram.
    Hist dndy("Rapidity distribution dn/dy of primary hadrons", 100,
	      -10., 10.);

    // Event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      // Reset event record and add q-qbar pair.
      event.reset();
      double mm = masslessQuarks ? 0 : pdt.m0(qid);
      double ee = smass / 2;
      double pp = sqrtpos(ee*ee - mm*mm);
      event.append(qid, 23, 101, 0, 0., 0., pp, ee, mm);
      event.append(-qid, 23, 0, 101, 0., 0., -pp, ee, mm);

      // Generate event.
      if (!pythia.next()) {
        cout << "Error: Event generation failed." << endl;
	break;
      }

      // Loop over particles.
      for (int i = 0; i < event.size(); ++i) {
        // Add primary hadron rapidities to histogram.
	int status = event[i].status();
	if (status > 80 && status < 90) {
          dndy.fill(event[i].y());
	}
      }
    }

    // Print statistics and histograms.
    pythia.stat();
    cout << dndy;

    // Add histogram to matplotlib output.
    hpl.add(dndy, "--," + plotColours[iMass], to_string_2dp(smass)
	    + " GeV string");
  }

  // Finalise and return
  hpl.plot();
  return 0;
}
