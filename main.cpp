#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"
#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

#include "tracker.h"
#include "TString.h"
#include <filesystem>

using Acts::UnitConstants::cm;
using std::numbers::pi;

int main(int argc, char *argv[]){
  auto trackingGeometryPtr = CreateTrackingGeometry();
  auto trackingGeometry = std::make_shared<Acts::TrackingGeometry>(*trackingGeometryPtr);

  // Random number generator config
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(ActsExamples::RandomNumbers::Config({42}));

  int nEvents = 1;

  ActsExamples::ParametricParticleGenerator::Config genCfg;
  genCfg.etaUniform = true;
  genCfg.phiMin = 0.;
  genCfg.phiMax = 0. + 0.0001;
  genCfg.thetaMin = pi + 1.0000e-1;
  genCfg.thetaMax = pi + 1.0001e-1;
  genCfg.phiMin = - pi;
  genCfg.phiMax = + pi;
  genCfg.thetaMin =  0 + 1.0000e-1;
  genCfg.thetaMax = pi - 1.0001e-1;

  // sec 3, positive x -> should be negative loc1 (local x)
  genCfg.phiMin = pi/2-pi/20;
  genCfg.phiMax = pi/2-pi/20+0.001;
  genCfg.thetaMin = pi/3;
  genCfg.thetaMax = pi/3 + 1.0001e-1;
  genCfg.pMin = 10.00;
  genCfg.pMax = 10.00+1e-10;
  genCfg.pTransverse = true;
  genCfg.pdg = Acts::ePionPlus;
  
  auto vertexGen = std::make_shared<ActsExamples::FixedPrimaryVertexPositionGenerator>();
  vertexGen->fixed = {0.0, 0.0, 0.0, 0.0};
  ActsExamples::EventGenerator::Generator gen{
      std::make_shared<ActsExamples::FixedMultiplicityGenerator>(1),
      vertexGen,
      std::make_shared<ActsExamples::ParametricParticleGenerator>(genCfg)};
  ActsExamples::EventGenerator::Config evgenCfg;
  evgenCfg.outputParticles = "particles";
  evgenCfg.outputVertices = "vertices";
  evgenCfg.generators = {gen};
  evgenCfg.randomNumbers = rnd;


  // Fatras config
  ActsExamples::FatrasSimulation::Config fatrasCfg;
  fatrasCfg.inputParticles = "particles";
  fatrasCfg.outputParticles = "final";
  fatrasCfg.outputSimHits = "simhits";
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.pMin = 0.05_GeV;
  fatrasCfg.magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 0.5*Acts::UnitConstants::T));
  fatrasCfg.randomNumbers = rnd;


  // Digitization config
  ActsExamples::DigitizationAlgorithm::Config digiCfg;
  digiCfg.inputSimHits = "simhits";
  digiCfg.randomNumbers = rnd;
  digiCfg.outputMeasurements = "measurements";
  digiCfg.surfaceByIdentifier = trackingGeometry->geoIdSurfaceMap();
  
  ActsExamples::DigiComponentsConfig digiConfig;
  digiConfig.smearingDigiConfig.push_back(ActsExamples::ParameterSmearingConfig{Acts::eBoundLoc0, ActsExamples::Digitization::Gauss(0.08)});
  digiConfig.smearingDigiConfig.push_back(ActsExamples::ParameterSmearingConfig{Acts::eBoundLoc1, ActsExamples::Digitization::Gauss(0.08)});
  std::vector<std::pair<Acts::GeometryIdentifier, ActsExamples::DigiComponentsConfig>> elements = { {Acts::GeometryIdentifier{}, digiConfig} };
  digiCfg.digitizationConfigs = Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>(elements);

  // Particle writer config
  ActsExamples::RootParticleWriter::Config particleWriterCfg;
  particleWriterCfg.inputParticles = "particles";
  particleWriterCfg.treeName = "particles";
  particleWriterCfg.filePath = "particles.root";

  // SimhitWriter config
  ActsExamples::RootSimHitWriter::Config simhitWriterCfg;
  simhitWriterCfg.inputSimHits = "simhits";
  simhitWriterCfg.filePath = "hits.root";

  ActsExamples::RootMeasurementWriter::Config measWriterCfg;
  measWriterCfg.inputMeasurements = "measurements";
  measWriterCfg.inputSimHits = "simhits";
  measWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  measWriterCfg.filePath = "measurements.root";
  measWriterCfg.surfaceByIdentifier = trackingGeometry->geoIdSurfaceMap();

  // Sequencer config
  ActsExamples::Sequencer::Config sequencerCfg;
  sequencerCfg.numThreads = 1;
  sequencerCfg.events = nEvents; 
  sequencerCfg.logLevel = Acts::Logging::ERROR;

  // Start sequencer
  ActsExamples::Sequencer sequencer(sequencerCfg);
  sequencer.addReader(std::make_shared<ActsExamples::EventGenerator>(evgenCfg, Acts::Logging::INFO));
  sequencer.addElement(std::make_shared<ActsExamples::FatrasSimulation>(fatrasCfg,Acts::Logging::INFO));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::DigitizationAlgorithm>(digiCfg, Acts::Logging::INFO));
  sequencer.addWriter(std::make_shared<ActsExamples::RootParticleWriter>(particleWriterCfg, Acts::Logging::INFO));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSimHitWriter>(simhitWriterCfg, Acts::Logging::INFO));
  sequencer.addWriter(std::make_shared<ActsExamples::RootMeasurementWriter>(measWriterCfg, Acts::Logging::INFO));


  sequencer.run();

  return 0;
}
