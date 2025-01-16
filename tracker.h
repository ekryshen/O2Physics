#ifndef tracker
#define tracker
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/NavigationLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Geant4/G4Material.hh"
#include "Geant4/G4NistManager.hh"

#include "TString.h"


//  native units: mm, GeV, c=1, e=1
using Acts::UnitConstants::cm;
using std::numbers::pi;

class MyDetectorElement : public Acts::DetectorElementBase {
public:
  MyDetectorElement(std::shared_ptr<const Acts::Transform3> transform, std::shared_ptr<Acts::Surface> surface, double thickness)
      : Acts::DetectorElementBase(), mTransform(transform), mSurface(surface), mThickness(thickness)
  {
  }
  const Acts::Transform3 &transform(const Acts::GeometryContext &gctx) const { return *mTransform; }
  const Acts::Surface &surface() const { return *mSurface; }
  Acts::Surface &surface() { return *mSurface; }
  double thickness() const { return mThickness; }

private:
  std::shared_ptr<const Acts::Transform3> mTransform = nullptr;
  std::shared_ptr<Acts::Surface> mSurface = nullptr;
  double mThickness = 0.;
};

std::vector<std::shared_ptr<MyDetectorElement>> detectorStore;
Acts::GeometryContext gctx;

Acts::TrackingGeometry* CreateTrackingGeometry(){
  // Create materials
  auto* nist = G4NistManager::Instance();
  G4Material* siMat = nist->FindOrBuildMaterial("G4_Si");
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");

  Acts::Geant4MaterialConverter converter;
  Acts::Material silicon = converter.material(*siMat);
  Acts::Material vacuum = converter.material(*worldMat);

  Acts::MaterialSlab matProp(silicon, 0.01*cm);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
  const auto volumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);

  int nLayers = 53;
  std::vector<int> padCount{20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31,
                            32, 32, 33, 33, 34, 35, 35, 36, 37, 38, 39, 40, 41, 41, 42, 43, 44, 45,
                            46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};

  Acts::TrackingVolumeVector vVolumes;
  std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3>> vVolumesBinning;

  for (int sec = 0; sec < 12 ; sec++) {
    Acts::LayerVector layVec;
    double phi0 = pi * (sec / 6. - 1);
    printf("phi0 = %f\n",phi0);
    for (int i = nLayers-1; i >=0; i--) {
      double y = 407. + ( i < 27 ? (i+0.5)*12. : (27*12. + (i - 27 +0.5)*18.));
      double x = padCount[i]*5.;
      double z = 1640.;
      // printf("%d %d %f\n", i, padCount[i], y);
      // const auto pBounds = std::make_shared<const Acts::RectangleBounds>(z, x); // half-x, half-y
      // Acts::Transform3 trafo = Acts::Transform3::Identity();
      // trafo.rotate(Eigen::AngleAxisd(pi/2., Acts::Vector3(0, 1, 0)));
      // trafo.rotate(Eigen::AngleAxisd(-phi0, Acts::Vector3(1, 0, 0)));
      // trafo.translate(Acts::Vector3(0., 0., y));

      const auto pBounds = std::make_shared<const Acts::RectangleBounds>(x, z); // half-x, half-y
      Acts::Transform3 trafo = Acts::Transform3::Identity();
      trafo.rotate(Eigen::AngleAxisd(pi/2., Acts::Vector3(1, 0, 0)));
      trafo.rotate(Eigen::AngleAxisd(pi/2., Acts::Vector3(0, 1, 0)));
      trafo.rotate(Eigen::AngleAxisd( phi0, Acts::Vector3(0, 1, 0)));
      trafo.translate(Acts::Vector3(0., 0., y));

      auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds);
      surface->assignSurfaceMaterial(std::move(surfaceMaterial));
      auto surArray = std::unique_ptr<Acts::SurfaceArray>(new Acts::SurfaceArray(surface));
      auto layer = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm);
      surface->associateLayer(*layer.get());
      layVec.push_back(layer);

      auto detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, 1._mm);
      surface->assignDetectorElement(*detElement.get());
      detectorStore.push_back(std::move(detElement));
    }

    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
    auto binningValue = ((sec>=2 && sec<=4) || (sec>=8 && sec<=10)) ? Acts::BinningValue::binY : Acts::BinningValue::binX;
    std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(gctx, layVec, -1300, 1300, Acts::BinningType::arbitrary, binningValue));

    auto boundsVol = std::make_shared<Acts::CylinderVolumeBounds>(200, 1400, 1700, pi/12., phi0); // rmin, rmax, halfz, half-phi, average-phi
    auto trackVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), boundsVol, volumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, Form("TPC%d",sec));
    vVolumes.push_back(trackVolume);
    vVolumesBinning.push_back(std::pair<Acts::TrackingVolumePtr, Acts::Vector3>(trackVolume, Acts::Vector3(cos(phi0),sin(phi0),0)));
  }

  // glue TPC sector volumes with each other
  for (int i=0; i<12; i++){
    auto mutableGlueVol1 = std::const_pointer_cast<Acts::TrackingVolume>(vVolumes[i]);
    auto mutableGlueVol2 = std::const_pointer_cast<Acts::TrackingVolume>(vVolumes[(i+1)%12]);
    mutableGlueVol2->glueTrackingVolume(gctx, Acts::BoundarySurfaceFace::index4, mutableGlueVol1.get(), Acts::BoundarySurfaceFace::index5);
  }

  // create pipe volume
  auto pipeBounds = std::make_shared<Acts::CylinderVolumeBounds>(0, 200., 1700.); // rmin, rmax, halfz
  auto pipeVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), pipeBounds, volumeMaterial, nullptr, nullptr, Acts::MutableTrackingVolumeVector{}, "PIPE");
  vVolumes.push_back(pipeVolume);

  // attach TPC volumes to the pipe outer surface
  auto pipeSurface = pipeVolume->boundarySurfaces().at(Acts::BoundarySurfaceFace::index2); // std::shared_ptr<const Acts::BoundarySurfaceT<Acts::TrackingVolume>>
  auto mutablePipeSurface = std::const_pointer_cast<Acts::BoundarySurfaceT<Acts::TrackingVolume>>(pipeSurface);
  auto binning = std::make_unique<const Acts::BinUtility>(12, -pi-pi/12, +pi-pi/12, Acts::closed, Acts::BinningValue::binPhi); // printf("%s\n",binning->toString().data());
  auto vArray = std::make_shared<const Acts::BinnedArrayXD<Acts::TrackingVolumePtr>>(vVolumesBinning, std::move(binning));
  mutablePipeSurface->attachVolumeArray(vArray, Acts::Direction::Positive);

  // create mother volume
  Acts::TrackingVolumeArrayCreator::Config volumeArrayCreatorCfg;
  Acts::TrackingVolumeArrayCreator volumeArrayCreator(volumeArrayCreatorCfg);
  auto volumeArray = volumeArrayCreator.trackingVolumeArray(gctx, vVolumes, Acts::BinningValue::binR);
  auto motherBounds = std::make_shared<Acts::CylinderVolumeBounds>(0, 1400., 1700.);
  auto motherVolume = std::make_shared<Acts::TrackingVolume>(Acts::Transform3::Identity(), motherBounds, volumeMaterial, nullptr, std::move(volumeArray), Acts::MutableTrackingVolumeVector{}, "DET");

  auto trackingGeometry = new Acts::TrackingGeometry(motherVolume);

  // Check geometry
  bool print_geometry_info = 1;
  if (print_geometry_info) {
    const Acts::TrackingVolume *highestTrackingVolume = trackingGeometry->highestTrackingVolume();
    printf("volumeId = %lu\n", highestTrackingVolume->geometryId().value());
    auto confinedVolumes = highestTrackingVolume->confinedVolumes();
    printf("confined volumes: %lu\n",confinedVolumes->arrayObjects().size());
    const Acts::LayerArray *confinedLayers = confinedVolumes->arrayObjects().at(12)->confinedLayers();
    for (const auto &layer : confinedLayers->arrayObjects())  {
      std::stringstream lid;
      lid << layer->geometryId();

      printf("  layerId = %s, thickness = %f type = %d\n", lid.str().data(), layer->thickness(), layer->layerType());
      if (layer->layerType()==-1) {
        const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
        // printf("   navigation layer%p\n", nlayer);
      } else {
        const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
        // printf("   plane %p\n", player);

      }
      if (!layer->surfaceArray())
        continue;
      for (const auto &surface : layer->surfaceArray()->surfaces()) {
        std::stringstream sid;
        sid << surface->geometryId();
        printf("       surfaceId = %s\n", sid.str().data());
      }
    } //for layers
  } // print_geometry_info


  return trackingGeometry;
}

#endif