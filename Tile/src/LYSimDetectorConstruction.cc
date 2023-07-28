#ifdef CMSSW_GIT_HASH
#include "MATHUSLA/Tile/interface/LYSimDetectorConstruction.hh"
#include "MATHUSLA/Tile/interface/LYSimDetectorMessenger.hh"
#include "MATHUSLA/Tile/interface/Materials.hh"
#include "MATHUSLA/Tile/interface/SurfaceProperty.hh"
#else
#include "LYSimDetectorConstruction.hh"
#include "LYSimDetectorMessenger.hh"
#include "Materials.hh"
#include "SurfaceProperty.hh"
#endif

#include <math.h>

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4Orb.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4EllipticalTube.hh"

using std::atan;
using std::cos;
using std::exp;
using std::sin;
using std::tan;
using namespace CLHEP;

LYSimPMTSD *LYSimDetectorConstruction::fPMTSD = NULL;

LYSimDetectorConstruction::LYSimDetectorConstruction()
    : G4VUserDetectorConstruction()
{
  fdetectorMessenger = new LYSimDetectorMessenger(this);

  ///////////////////////////////////////////////////////////////////////////////
  // dimensions
  ///////////////////////////////////////////////////////////////////////////////

  _tilex = 10 * mm;
  _tiley = 10 * mm;
  _tilez = 50 * mm; 
  _tile_x1 = 0.0 * mm; //extrusion1 off center
  _tile_x2 = 0.0 * mm; //extrusion2 off center
  wrapgap = 0.1 * mm; // gap between wrap and scintillator, disabled if using cladding wrap
  wrapthickness = 0.1 * mm; // thickness of wrap

  // SiPM property
  // _sipm_deadwidth = 0.2 * mm;
  // _sipm_x = 1.4 * mm;
  // _sipm_y = 1.4 * mm;
  _sipm_z = 0.4 * mm;
  // _sipm_rimwidth = 0.3 * mm;
  // _sipm_glasswidth = 0.1 * mm;
  // _sipm_standz = 0.3 * mm;
  // _sipm_x = 10 * mm;
  // _sipm_y = 10 * mm;
  // Default Dimple settings
  // _dimple_type = SPHERICAL;
  // _dimple_indent = 1.0 * mm;
  // _dimple_radius = 6.0 * mm; // 3.4409*mm
  // Default Hole settings
  // _pcb_radius = 2.5;

  // Fiber/extrusion related
  _cladlayer  = 1; //number of fiber clad layers 
  _holeshape = 0; //0 circle; 1 square; 2 el
  _WLSfiberR = 0.7*mm;
  _WLSfiber_clad_thick = 0.05*mm;
  _WLSfiber_clad2_thick = 0.05*mm;
  _hole_radius = 1*mm; // hole raidus
  _hole_x1 = 0*mm; // hole position 
  _hole_x2 = _hole_x1; // a second hole position
  _WLSfiberZ = 3.2*m; // length of fiber
  _WLS_zoff = 1.7*m; // fiber offset
  _WLSfiberZ = _tilez*2.5;
  _WLS_zoff = 0;
  _WLS_xoff = 0;
  // new setting
  _holeshape = 2; // 0 circle; 1 square; 2 el
  _hole_radius = (_WLSfiberR + _WLSfiber_clad_thick + _WLSfiber_clad2_thick) * 1.01; // 1.0*mm;
  _hole_radius = 3 * mm;
  _hole_x1 = 0 * mm;
  _hole_x2 = 13 * mm;
  _WLSfiberZ = _tilez * 1.5;

  ///////////////////////////////////////////////////////////////////////////////
  // material
  ///////////////////////////////////////////////////////////////////////////////
  _absmult = 800;  // scaling of absorption, unit is mm
  _ScintiN = 100;  // scintillation scaling, unit / MeV
  _tiledecay = 5;  // ns // tuning to 11.53 ns overall
  _tilerise = 0.5; // ns
  _tile_alpha = 0.01;
  _dimple_alpha = 0.1;

  mfiber  = Make_Y11();
  mfiber_clad = Make_acrylic(); //cald 1
  mfiber_clad2 = Make_Fluor_acrylic(); //clad 2
  _fiber_decaytime = 11.5; //ns
  fcoating = Make_Coating();  // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.

  fBialkali = Make_Bialkali(); // SiPM material
  fEpoxy = Make_Epoxy(); // wrapper material
  fAir = Make_Custom_Air(); // well, air 
  fEJ200 = Make_EJ200(); // Scintillator material
  fResin = Make_Resin(); // Fill the gap 
  SetTileAbsMult(_absmult);
  SetTileScintillation(_ScintiN); // N/keV
  SetTileDecaytime(_tiledecay);
  SetTileRisetime(_tilerise);
  SetFiberdecaytime( _fiber_decaytime );


  ///////////////////////////////////////////////////////////////////////////////
  // surface
  ///////////////////////////////////////////////////////////////////////////////
  _handwrap   = false; // if _handwrap==True, no air gap between wrap and scintillator 
  _coverends  = true;
  // wrapper
  fESROpSurface = MakeS_RoughMirror(); // surface of the side wrap 
  fTiO2Surface = MakeS_TiO2Surface(); // Extrusion (TiO2 Coating) Surface
  ffrontendSurface = MakeS_VAbsorbing();

  // _sipm_eff = 1;
  // Defining surface list.
  // fIdealPolishedOpSurface = MakeS_IdealPolished();
  fTileBulkSurface = MakeS_RoughInterface(_tile_alpha);
  fSiPMSurface3 = MakeS_SiPM();
  fSiPMSurface4 = MakeS_SiPM();
  // _claddirt = 0.0;
  // opSurface =  MakeS_Polished();
  //opSurface = MakeS_IdealPolished();
  _wrap_reflect = 0.985;
  _frontend_reflect = 0.985;
  SetWrapReflect(_wrap_reflect);
  SetfrontendReflect(_frontend_reflect);
}

void LYSimDetectorConstruction::UpdateGeometry()
{
  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();

  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

LYSimDetectorConstruction::~LYSimDetectorConstruction()
{
  if (fdetectorMessenger)
  {
    delete fdetectorMessenger;
  }
}

G4VPhysicalVolume *
LYSimDetectorConstruction::Construct()
{

  static const bool checkOverlaps = true;
  double verbose=11;
  ///////////////////////////////////////////////////////////////////////////////
  // World volume
  ///////////////////////////////////////////////////////////////////////////////
  G4Box *solidWorld = new G4Box("World", WorldHalfX(), WorldHalfY(), WorldHalfZ());

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");

  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWorld, "World", 0, false, 0, checkOverlaps);

  ///////////////////////////////////////////////////////////////////////////////
  // wrapping or clad
  ///////////////////////////////////////////////////////////////////////////////
  G4VSolid *solidHoleBound;
  // Set hole shape
  if (_holeshape == 0)
    solidHoleBound = new G4Tubs("TileHole", 0, _hole_radius, _tilez, 0, 2 * pi);
  else if (_holeshape == 1)
    solidHoleBound = new G4Box("TileHole", _hole_radius, _hole_radius, _tilez);
  else if (_holeshape == 2)
    solidHoleBound = new G4EllipticalTube("TileHole", 1.2, _hole_radius, _tilez);

  G4LogicalVolume *logicWrap;
  G4LogicalVolume *logicWrapface;
  G4VPhysicalVolume *physWrap3;
  G4VPhysicalVolume *physWrap4;
  // set wrap at side, no end covering
  if (_handwrap)
  {
    // handwrap, with air gap
    G4VSolid *solidWrap0 = ConstructHollowWrapSolid();
    G4VSolid *solidWrap1 = new G4SubtractionSolid("solidWrap", solidWrap0, solidHoleBound, 0, G4ThreeVector(_hole_x1, 0, 0));
    // second hole
    // G4VSolid *solidWrap = new G4SubtractionSolid("solidWrap", solidWrap1, solidHoleBound, 0, G4ThreeVector(_hole_x2, 0, 0));
    // logicWrap = new G4LogicalVolume(solidWrap, fEpoxy, "Wrap");
    logicWrap = new G4LogicalVolume(solidWrap1, fEpoxy, "Wrap"); // comment this line when second hole enabled
    G4VPhysicalVolume *physWrap = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWrap, "Wrap", logicWorld, false, 0, checkOverlaps);
    G4LogicalSkinSurface *WrapSurface = new G4LogicalSkinSurface("WrapSurface", logicWrap, fESROpSurface);
  }
  else
  {
    // clad, no air gap
    G4VSolid *solidExtrusion = ConstructHollowWrapCladSolid();
    logicWrap = new G4LogicalVolume(solidExtrusion, fcoating, "Extrusion");
    G4VPhysicalVolume *physWrap = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicWrap, "Wrap", logicWorld, false, 0, checkOverlaps);
    G4LogicalSkinSurface *WrapSurface = new G4LogicalSkinSurface("WrapSurface", logicWrap, fTiO2Surface);//surface use clading 
  }
  
  ///////////////////////////////////////////////////////////////////////////////
  // extruded scintillator
  ///////////////////////////////////////////////////////////////////////////////
  G4VSolid *solidTile = ConstructTrapazoidSolid("TileTrap", _tilex, _tiley, _tilez, _tile_x1, _tile_x2);
  G4VSolid *tileBulk1 = new G4SubtractionSolid("TileSolid_Bulk", solidTile, solidHoleBound, 0, G4ThreeVector(_hole_x1, 0, 0));
  // second hole ----------
  // G4VSolid *tileBulk = new G4SubtractionSolid("TileSolid_Bulk", tileBulk1, solidHoleBound, 0, G4ThreeVector(_hole_x2, 0, 0));
  // G4LogicalVolume *logicTileBulk = new G4LogicalVolume(tileBulk, fEJ200, "TileBulkLogic");
  //----------------
  G4LogicalVolume *logicTileBulk = new G4LogicalVolume(tileBulk1, fEJ200, "TileBulkLogic"); // comment this line when second hole enabled
  G4VPhysicalVolume *physTileBulk = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicTileBulk, "TileBulkPhysic", logicWorld, false, 0, checkOverlaps);

  // set wrap at front/back ends
  if(_coverends){

    G4VSolid *wrapface = ConstructTrapazoidSolid("Wrapface", _tilex, _tiley, wrapthickness, 0, 0);
    G4VSolid *solidWrapface = new G4SubtractionSolid("solidWrapface", wrapface, solidHoleBound, 0, G4ThreeVector(_hole_x1, 0, 0));
    logicWrapface = new G4LogicalVolume(solidWrapface, fEpoxy, "Wrapface");
    G4LogicalSkinSurface *WrapfaceSurface = new G4LogicalSkinSurface("WrapfaceSurface", logicWrapface, ffrontendSurface); // define surface at end 

    physWrap3 = new G4PVPlacement(0, G4ThreeVector(0, 0, _tilez * 0.5 + wrapthickness * 0.5), logicWrapface, "Wrapface3", logicWorld, false, 0, checkOverlaps);
    physWrap4 = new G4PVPlacement(0, G4ThreeVector(0, 0, -_tilez * 0.5 - wrapthickness * 0.5), logicWrapface, "Wrapface4", logicWorld, false, 0, checkOverlaps);
    
    new G4LogicalBorderSurface("surfacewrapout1", physWrap3, physTileBulk, MakeS_Absorbing());
    new G4LogicalBorderSurface("surfacewrapout2", physTileBulk, physWrap3, MakeS_Absorbing());
    new G4LogicalBorderSurface("surfacewrapout3", physTileBulk, physWrap4, MakeS_Absorbing());
    new G4LogicalBorderSurface("surfacewrapout4", physWrap4, physTileBulk, MakeS_Absorbing());
    
  }


  ///////////////////////////////////////////////////////////////////////////////
  // fiber. single layer clad
  ///////////////////////////////////////////////////////////////////////////////
  assert(_hole_radius >= _WLSfiberR + _WLSfiber_clad_thick + _WLSfiber_clad2_thick);
  G4VSolid *solidWLSfiber = new G4Tubs("WLSFiber", 0., _WLSfiberR, _WLSfiberZ * 0.5, 0., 2 * pi);
  G4VSolid *solidWLSfiber_clad = new G4Tubs("WLSFiber_clad", _WLSfiberR, _WLSfiberR + _WLSfiber_clad_thick, _WLSfiberZ * 0.5, 0., 2 * pi);
  G4LogicalVolume *logicWLSfiber = new G4LogicalVolume(solidWLSfiber, mfiber, "logicWLSfiber");
  G4LogicalVolume *logicWLSfiber_clad = new G4LogicalVolume(solidWLSfiber_clad, mfiber_clad, "logicWLSfiber_clad");

  G4VPhysicalVolume *physWLSfiber = new G4PVPlacement(0, G4ThreeVector(_hole_x1 + _WLS_xoff, 0, -_WLS_zoff), logicWLSfiber, "PhyhWLSfiber_core", logicWorld, false, 0, checkOverlaps);

  G4VPhysicalVolume *physWLSfiber_clad = new G4PVPlacement(0, G4ThreeVector(_hole_x1 + _WLS_xoff, 0, -_WLS_zoff), logicWLSfiber_clad, "PhyhWLSfiber_cald", logicWorld, false, 0, checkOverlaps);

  ///////////////////////////////////////////////////////////////////////////////
  // fiber. second layer clad
  ///////////////////////////////////////////////////////////////////////////////
  G4LogicalVolume *logicWLSfiber_clad2;
  G4VSolid *solidWLSfiber_clad2;
  std::cout << "============  add a second cladding! =================" << std::endl;
  solidWLSfiber_clad2 = new G4Tubs("WLSFiber_clad2", _WLSfiberR + _WLSfiber_clad_thick, _WLSfiberR + _WLSfiber_clad_thick + _WLSfiber_clad2_thick, _WLSfiberZ * 0.5, 0., 2 * pi);
  logicWLSfiber_clad2 = new G4LogicalVolume(solidWLSfiber_clad2, mfiber_clad2, "logicWLSfiber_clad2");

  if (_cladlayer == 2)
  {
    G4VPhysicalVolume *physWLSfiber_clad2 = new G4PVPlacement(0, G4ThreeVector(_hole_x1 + _WLS_xoff, 0, -_WLS_zoff), logicWLSfiber_clad2, "PhyhWLSfiber_cald2", logicWorld, false, 0, checkOverlaps);
  }

  ///////////////////////////////////////////////////////////////////////////////
  // gap between fiber and hole
  ///////////////////////////////////////////////////////////////////////////////
  double holeinner = _WLSfiberR + _WLSfiber_clad_thick;
  if (_cladlayer == 2)
    holeinner += _WLSfiber_clad2_thick;
  // circle
  G4VSolid *solidHole;
  if (_holeshape == 0)
    solidHole = new G4Tubs("TileHole", 0, _hole_radius, _tilez * 0.5, 0, 2 * pi);
  else if (_holeshape == 1)
    solidHole = new G4Box("TileHole", _hole_radius, _hole_radius, _tilez * 0.5);
  else if (_holeshape == 2)
    solidHole = new G4EllipticalTube("TileHole",  1.2, _hole_radius, _tilez * 0.5);

  G4VSolid *solidWLSfiberFrame = new G4Tubs("WLSFiberFrame", 0., holeinner, _WLSfiberZ * 0.5, 0., 2 * pi);
  G4VSolid *solidHole_subs = new G4SubtractionSolid("TileHole_Subs", solidHole, solidWLSfiberFrame, 0, G4ThreeVector(_WLS_xoff, 0, 0));
  //G4LogicalVolume *logicHole = new G4LogicalVolume(solidHole_subs, fAir, "Hole"); // normal air gap
  G4LogicalVolume *logicHole = new G4LogicalVolume(solidHole_subs, fResin, "Hole");//testing gap material 

  G4VPhysicalVolume *physHole = new G4PVPlacement(0, G4ThreeVector(_hole_x1, 0, 0), logicHole, "Hole", logicWorld, false, 0, checkOverlaps);

  /*
    if(_cladlayer==2){
          new G4LogicalBorderSurface("surfacefiberOut", physWLSfiber,             physWLSfiber_clad, opSurface);
          new G4LogicalBorderSurface("surfacefiberin", physWLSfiber_clad,             physWLSfiber, opSurface);
          new G4LogicalBorderSurface("surfaceclad1Out", physWLSfiber_clad,             physWLSfiber_clad2, opSurface);
          new G4LogicalBorderSurface("surfaceclad1in", physWLSfiber_clad2,             physWLSfiber_clad, opSurface);
          new G4LogicalBorderSurface("surfaceClad2Out", physWLSfiber_clad2,             physHole, opSurface);
          new G4LogicalBorderSurface("surfaceClad2In", physHole, physWLSfiber_clad2,             opSurface);
    }else{
          new G4LogicalBorderSurface("surfacefiber1Out", physWLSfiber,             physWLSfiber_clad, opSurface);
          new G4LogicalBorderSurface("surfacefiber1in", physWLSfiber_clad,             physWLSfiber, opSurface);
          new G4LogicalBorderSurface("surfaceClad1Out", physWLSfiber_clad,             physHole, opSurface);
          new G4LogicalBorderSurface("surfaceClad1In",  physHole, physWLSfiber_clad,             opSurface);

    }
  */
/* //remove clad
  G4OpticalSurface *surface_IdealnonPolished = new G4OpticalSurface("IdealnonPolished", glisur, ground, dielectric_dielectric, 1);
  double surface_IdealnonPolishedphoE2[2] = {1.0 * eV, 8.0 * eV};
  double surface_IdealnonPolishedreflectivity[2] = {1 - _claddirt, 1 - _claddirt};
  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", surface_IdealnonPolishedphoE2, surface_IdealnonPolishedreflectivity, 2);
  // table->AddProperty( "TRANSMITTANCE", surface_IdealnonPolishedphoE2, surface_IdealnonPolishedreflectivity, 2 );
  surface_IdealnonPolished->SetMaterialPropertiesTable(table);

  new G4LogicalSkinSurface("fibercore", logicWLSfiber, opSurface);
  // new G4LogicalSkinSurface( "fiberclad", logicWLSfiber_clad, surface_IdealnonPolished );

  new G4LogicalBorderSurface("surfaceClad1Out", physWLSfiber_clad,     physHole, surface_IdealnonPolished);
  new G4LogicalBorderSurface("surfaceClad1In", physHole, physWLSfiber_clad,     surface_IdealnonPolished);
  // end here
*///------------remove clad
  // G4OpticalSurface* surface_p =MakeS_TiO2Surface();
  // new G4LogicalBorderSurface("surfaceCladin", physWLSfiber_clad, physHole, surface_p);
  /*new G4LogicalBorderSurface("surfacecladout", physWLSfiber_clad, physHole, opSurface);
  new G4LogicalBorderSurface("surfacecladin", physHole, physWLSfiber_clad, opSurface);
  new G4LogicalBorderSurface("surfacecoreout", physWLSfiber_clad, physWLSfiber, opSurface);
  new G4LogicalBorderSurface("surfacecorein", physWLSfiber, physWLSfiber_clad, opSurface);
  */
  ///////////////////////////////////////////////////////////////////////////////
  // realistic SiPM
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  // Simple version of SiPM
  ///////////////////////////////////////////////////////////////////////////////

  const G4ThreeVector SiPMOffset_chan3(_hole_x1 + _WLS_xoff, 0, _WLSfiberZ * 0.5 - _WLS_zoff + 0.5 * _sipm_z + 0.01 * mm);
  const G4ThreeVector SiPMOffset_chan4(_hole_x1 + _WLS_xoff, 0, -_WLSfiberZ * 0.5 - _WLS_zoff - 0.5 * _sipm_z - 0.01 * mm); // air gap

  // G4Tubs* solidSiPMInnerBox = new G4Tubs( "solidSiPMInnerBox", 0., _WLSfiberR, _sipm_z*0.5, 0., 2*pi);
  G4VSolid *solidSiPMInnerBox = new G4Box("solidSiPMInnerBox", _WLSfiberR, _WLSfiberR, _sipm_z * 0.5);

  G4LogicalVolume *logicSiPM = new G4LogicalVolume(solidSiPMInnerBox, fBialkali, "SiPM");

  G4VPhysicalVolume *physSiPM_chan3 = new G4PVPlacement(0, SiPMOffset_chan3, logicSiPM, "physSiPM_chan3", logicWorld, false, 0, checkOverlaps);
  G4VPhysicalVolume *physSiPM_chan4 = new G4PVPlacement(0, SiPMOffset_chan4, logicSiPM, "physSiPM_chan4", logicWorld, false, 0, checkOverlaps);

  new G4LogicalBorderSurface("SiPMSurface3_out", physWorld, physSiPM_chan3, fSiPMSurface3);
  new G4LogicalBorderSurface("SiPMSurface4_out", physWorld, physSiPM_chan4, fSiPMSurface4);
  new G4LogicalBorderSurface("SiPMSurface3_in", physSiPM_chan3, physWorld, fSiPMSurface3);
  new G4LogicalBorderSurface("SiPMSurface4_in", physSiPM_chan4, physWorld, fSiPMSurface4);

  //    if(_cladlayer==2){
  //        new G4LogicalBorderSurface("SiPMSurface3_out1", physWLSfiber_clad, physSiPM_chan3, fSiPMSurface3);
  //        new G4LogicalBorderSurface("SiPMSurface4_out1", physWLSfiber_clad, physSiPM_chan4, fSiPMSurface4);
  //        new G4LogicalBorderSurface("SiPMSurface3_out2", physWLSfiber_clad2, physSiPM_chan3, fSiPMSurface3);
  //        new G4LogicalBorderSurface("SiPMSurface4_out2", physWLSfiber_clad2, physSiPM_chan4, fSiPMSurface4);
  //    }else{
  //        new G4LogicalBorderSurface("SiPMSurface3_out1", physWLSfiber_clad, physSiPM_chan3, fSiPMSurface3);
  //        new G4LogicalBorderSurface("SiPMSurface4_out1", physWLSfiber_clad, physSiPM_chan4, fSiPMSurface4);
  //    }

  G4OpticalSurface *FiberEndSurface = new G4OpticalSurface("FiberEndSurface",                                   unified,                                   ground,                                   dielectric_dielectric,                                   0.02);
  new G4LogicalBorderSurface("FiberEndSurface1", physWLSfiber, physWorld, FiberEndSurface);
  new G4LogicalBorderSurface("FiberEndSurface2", physWorld, physWLSfiber, FiberEndSurface);
  new G4LogicalBorderSurface("FiberEndSurface3", physWLSfiber_clad, physWorld, FiberEndSurface);
  new G4LogicalBorderSurface("FiberEndSurface4", physWorld, physWLSfiber_clad, FiberEndSurface);

  ///////////////////////////////////////////////////////////////////////////////
  // Defining surfaces
  ///////////////////////////////////////////////////////////////////////////////

  // Tile surfaces
  G4LogicalBorderSurface *TileBulkSurface =
      new G4LogicalBorderSurface("TileBulkSurface", physTileBulk, physWorld, fTileBulkSurface);

  // Other optical surfaces

  // Setting the sensitive detector
  if (!fPMTSD)
  {
    fPMTSD = new LYSimPMTSD("/LYSimPMT");
    G4SDManager *sdman = G4SDManager::GetSDMpointer();
    sdman->AddNewDetector(fPMTSD);
  }
  logicSiPM->SetSensitiveDetector(fPMTSD);

  G4Colour white(1.00, 1.00, 1.00);   // white
  G4Colour gray(0.50, 0.50, 0.50);    // gray
  G4Colour black(0.00, 0.00, 0.00);   // black
  G4Colour red(1.00, 0.00, 0.00);     // red
  G4Colour green(0.00, 1.00, 0.00);   // green
  G4Colour blue(0.00, 0.00, 1.00);    // blue
  G4Colour cyan(0.00, 1.00, 1.00);    // cyan
  G4Colour air(0.90, 0.94, 1.00);     // cyan
  G4Colour magenta(1.00, 0.00, 1.00); // magenta
  G4Colour yellow(1.00, 1.00, 0.00);  // yellow
  G4Colour brass(0.80, 0.60, 0.40);   // brass
  G4Colour brown(0.70, 0.40, 0.10);   // brown

  // Visual attributes
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *SiPMVisAtt = new G4VisAttributes(red);
  SiPMVisAtt->SetForceSolid(true);
  SiPMVisAtt->SetVisibility(true);
  // logicSiPM->SetVisAttributes( SiPMVisAtt );
  logicSiPM->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *ResinVisAtt = new G4VisAttributes(cyan);
  ResinVisAtt->SetForceWireframe(true);
  ResinVisAtt->SetVisibility(true);
  logicHole->SetVisAttributes(ResinVisAtt);

  G4VisAttributes *TileVisAtt = new G4VisAttributes(G4Colour(0, 0, 1.));
  TileVisAtt->SetForceWireframe(true);
  TileVisAtt->SetVisibility(true);
  //logicTileBulk->SetVisAttributes(G4VisAttributes::Invisible);
  logicTileBulk->SetVisAttributes(TileVisAtt);

  G4VisAttributes *fiberVisAtt = new G4VisAttributes(G4Colour(1, 0, 0));
  fiberVisAtt->SetForceWireframe(true);
  fiberVisAtt->SetVisibility(true);
  fiberVisAtt->SetForceAuxEdgeVisible(true);
  logicWLSfiber->SetVisAttributes(fiberVisAtt);

  G4VisAttributes *fiberVisAtt2 = new G4VisAttributes(G4Colour(0.5, 0.5, 0));
  fiberVisAtt2->SetForceWireframe(true);
  fiberVisAtt2->SetVisibility(true);
  fiberVisAtt2->SetForceAuxEdgeVisible(true);
  logicWLSfiber_clad->SetVisAttributes(fiberVisAtt2);
  if (_cladlayer == 2)
    logicWLSfiber_clad2->SetVisAttributes(fiberVisAtt2);

  G4VisAttributes *WrapVisAtt = new G4VisAttributes(gray);
  WrapVisAtt->SetForceWireframe(true);
  WrapVisAtt->SetVisibility(true);
  logicWrap->SetVisAttributes(WrapVisAtt);

  return physWorld;
}

G4VSolid *
LYSimDetectorConstruction::ConstructTrapazoidSolid(
    const G4String &name,
    double x,
    double y,
    double z,
    double indent_x1,
    double indent_x2) const
{
  const G4ThreeVector corners[8] = {
      G4ThreeVector(-x / 2, -y / 2, -z / 2),
      G4ThreeVector(x / 2, -y / 2, -z / 2),
      G4ThreeVector(-x / 2 + indent_x2, y / 2, -z / 2),
      G4ThreeVector(x / 2 - indent_x1, y / 2, -z / 2),
      G4ThreeVector(-x / 2, -y / 2, z / 2),
      G4ThreeVector(x / 2, -y / 2, z / 2),
      G4ThreeVector(-x / 2 + indent_x2, y / 2, z / 2),
      G4ThreeVector(x / 2 - indent_x1, y / 2, z / 2)};

  return new G4Trap(name, corners);
}

G4VSolid *
LYSimDetectorConstruction::ConstructHollowWrapSolid() const
{
  G4VSolid *wrapOuter = ConstructTrapazoidSolid("WrapOuter", _tilex + 2 * wrapgap + 2 * wrapthickness, _tiley + 2 * wrapgap + 2 * wrapthickness, _tilez + 2 * wrapgap + 2 * wrapthickness, 0, 0);
  G4VSolid *wrapInner = ConstructTrapazoidSolid("WrapInner", _tilex + 2 * wrapgap, _tiley + 2 * wrapgap, _tilez * 1.5, 0, 0);
  G4VSolid *wrapbox = new G4SubtractionSolid("WrapBox", wrapOuter, wrapInner);

  return wrapbox;
}

G4VSolid *
LYSimDetectorConstruction::ConstructHollowWrapCladSolid() const
{
  G4VSolid *wrapOuter = ConstructTrapazoidSolid("WrapOuter", _tilex + 2 * wrapthickness, _tiley + 2 * wrapthickness, _tilez, 0, 0);
  G4VSolid *wrapInner = ConstructTrapazoidSolid("WrapInner", _tilex, _tiley, _tilez * 1.5, 0, 0);
  G4VSolid *wrapbox = new G4SubtractionSolid("WrapBox", wrapOuter, wrapInner);
  return wrapbox;
}

G4VSolid *
LYSimDetectorConstruction::ConstructDimpleSubtract() const
{
  switch (_dimple_type)
  {
  case SPHERICAL:
    return ConstructSphereDimpleSolid();
  case FLAT_DOME:
    return ConstructFlatDomeDimpleSolid();
  case CYLINDRICAL:
    return ConstructCylindricalDimpleSolid();
  case ELLIPSOID:
    return ConstructEllipsoidDimpleSolid();
  default:
    return ConstructSphereDimpleSolid();
  }
}

G4VSolid *
LYSimDetectorConstruction::ConstructSphereDimpleSolid() const
{
  return new G4Orb("DimpleSolid", GetDimpleSphereRadius());
}

G4VSolid *
LYSimDetectorConstruction::ConstructFlatDomeDimpleSolid() const
{
  G4VSolid *torus = new G4Torus("DimpleTorus", 0, _dimple_indent, _dimple_radius - _dimple_indent, 0, 2 * pi);
  G4VSolid *tub = new G4Tubs("DimpleTub", 0, _dimple_radius - _dimple_indent, _dimple_indent, 0, 2 * pi);

  return new G4UnionSolid("DimpleSolid", torus, tub);
}

G4VSolid *
LYSimDetectorConstruction::ConstructCylindricalDimpleSolid() const
{
  return new G4Tubs("DimpleSolid", 0.0, _dimple_radius, _dimple_indent, 0, 2 * pi);
}

G4VSolid *
LYSimDetectorConstruction::ConstructEllipsoidDimpleSolid() const
{
  return new G4Ellipsoid("DimpleSolid", _dimple_radius, _dimple_radius, _dimple_indent);
}

G4ThreeVector
LYSimDetectorConstruction::CalcDimpleOffset() const
{
  switch (_dimple_type)
  {
  case FLAT_DOME:
    return G4ThreeVector(0, 0, 0.5 * _tilez);
  case CYLINDRICAL:
    return G4ThreeVector(0, 0, 0.5 * _tilez);
  case ELLIPSOID:
    return G4ThreeVector(0, 0, 0.5 * _tilez);
  default:
    return G4ThreeVector(0, 0, 0.5 * _tilez + GetDimpleSphereRadius() - _dimple_indent);
  }
}

// Additional geometry factors
double
LYSimDetectorConstruction::WorldHalfX() const
{
  return 1 * m;
}

double
LYSimDetectorConstruction::WorldHalfY() const
{
  return 1 * m;
}

double
LYSimDetectorConstruction::WorldHalfZ() const
{
  return 40 * m;
}

double
LYSimDetectorConstruction::LocalTileZ(const double x, const double y) const
{
  if (x < _tile_x1 + _hole_radius && x > _tile_x1 - _hole_radius)
  {
    return _tilez - 2 * sqrt(pow(_hole_radius, 2) - pow(x - _tile_x1, 2));
  }
  else if (x < _tile_x2 + _hole_radius && x > _tile_x2 - _hole_radius)
  {
    return _tilez - 2 * sqrt(pow(_hole_radius, 2) - pow(x - _tile_x2, 2));
  }
  else
  {
    return _tilez;
  }
}

double
LYSimDetectorConstruction::LocalTileZFlatDome(const double r) const
{
  if (r < _dimple_radius - _dimple_indent)
  {
    return _tilez - _dimple_indent;
  }
  else
  {
    const double dif = r - _dimple_radius + _dimple_indent;
    return _tilez - sqrt(_dimple_indent * _dimple_indent - dif * dif);
  }
}

double
LYSimDetectorConstruction::LocalTileZCylindrical(const double r) const
{
  return _tilez - _dimple_indent;
}

double
LYSimDetectorConstruction::LocalTileZEllipsoid(const double r) const
{
  return _tilez - _dimple_indent * sqrt(
                                       1 - r * r / _dimple_radius / _dimple_radius);
}

double
LYSimDetectorConstruction::LocalTileZSpherical(const double r) const
{
  const double big_r = GetDimpleSphereRadius();
  return _tilez - (sqrt(big_r * big_r - r * r) - (big_r - _dimple_indent));
}

double
LYSimDetectorConstruction::GetDimpleSphereRadius() const
{
  return 0.5 * ((_dimple_radius * _dimple_radius) / _dimple_indent + _dimple_indent);
}

void LYSimDetectorConstruction::SetTileAbsMult(const double mult)
{
  _absmult = mult;
  std::cout << "========== SetTileAbsMult to " << mult * mm << " mm ===================" << std::endl;
  Set_Flat_AbsLength(fEJ200, _absmult);
}

void LYSimDetectorConstruction::SetTileScintillation(const double mult)
{
  _ScintiN = mult;
  std::cout << "========== Update_EJ200_Scinti " << _ScintiN << " / MeV ===================" << std::endl;
  Update_EJ200_Scinti(fEJ200, _ScintiN);
}

void LYSimDetectorConstruction::SetTileDecaytime(const double mult)
{
  _tiledecay = mult;
  std::cout << "========== Set TileDecaytime " << _tiledecay << " ns ===================" << std::endl;
  Update_EJ200_FASTTIMECONSTANT(fEJ200, _tiledecay);
}

void LYSimDetectorConstruction::SetTileRisetime(const double mult)
{
  _tilerise = mult;
  std::cout << "========== Set TileRisetime " << _tilerise << " ns ===================" << std::endl;
  Update_EJ200_RISETIME(fEJ200, _tilerise);
}

void LYSimDetectorConstruction::SetFiberdecaytime(const double mult)
{
  _fiber_decaytime = mult;
  Update_Fiber_timeconstant(mfiber, mult);
  std::cout << "========== Update_Fiber_timeconstant to " << mult * ns << " ns ===================" << std::endl;
}

void LYSimDetectorConstruction::SetFiberattenu(const double mult)
{
  Set_Flat_AbsLength(mfiber, mult);
  std::cout << "========== SetFiberattenu to " << mult * mm << " mm ===================" << std::endl;
  // std::cout<<"========== Fiber abs length is perfect, don't update ==================="<<std::endl;
}

void LYSimDetectorConstruction::SetGaprefrac_index(const double mult)
{
  std::cout << "========== SetGaprefrac_index to " << mult << " ===================" << std::endl;
  Update_refrac_index(fResin, mult);
}

void LYSimDetectorConstruction::SetClad_refrac_index(const double mult1, const double mult2)
{
  std::cout << "========== SetClad_refrac_index to " << mult1 << " " << mult2 << " ===================" << std::endl;
  Update_refrac_index(mfiber_clad, mult1);
  Update_refrac_index(mfiber_clad2, mult2);
}

void LYSimDetectorConstruction::SetWrapReflect(const double r)
{
  // Add entries into properties table
  _wrap_reflect = r;
  static const unsigned nentries = 2;
  static double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  double reflectivity[nentries] = {_wrap_reflect, _wrap_reflect};

  G4MaterialPropertiesTable *table = fESROpSurface->GetMaterialPropertiesTable();
  if (table)
  {
    table->RemoveProperty("REFLECTIVITY");
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  }
  else
  {
    table = new G4MaterialPropertiesTable();
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
    fESROpSurface->SetMaterialPropertiesTable(table);
  }

  G4MaterialPropertiesTable* table2 = fTiO2Surface->GetMaterialPropertiesTable();
  if( table2 ){
    table2->RemoveProperty( "REFLECTIVITY" );
    table2->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table2 = new G4MaterialPropertiesTable();
    table2->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fTiO2Surface->SetMaterialPropertiesTable( table2 );
  }
}

void LYSimDetectorConstruction::SetfrontendReflect(const double r)
{
  // Add entries into properties table
  _frontend_reflect = r;
  static const unsigned nentries = 2;
  static double phoE[nentries] = {1.0 * eV, 6.0 * eV};
  double reflectivity[nentries] = {_frontend_reflect, _frontend_reflect};
  G4MaterialPropertiesTable *table = ffrontendSurface->GetMaterialPropertiesTable();
  if (table)
  {
    table->RemoveProperty("REFLECTIVITY");
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  }
  else
  {
    table = new G4MaterialPropertiesTable();
    table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
    ffrontendSurface->SetMaterialPropertiesTable(table);
  }

}

void LYSimDetectorConstruction::SetSiPMReflect(const double r)
{
  // Add entries into properties table
  _sipm_eff = r;
  G4double p_mppc[62] = {1.38 * eV, 1.39 * eV, 1.43 * eV, 1.43 * eV, 1.47 * eV, 1.48 * eV, 1.51 * eV, 1.53 * eV, 1.55 * eV, 1.59 * eV, 1.64 * eV, 1.68 * eV, 1.69 * eV, 1.72 * eV, 1.77 * eV, 1.78 * eV, 1.87 * eV, 1.92 * eV, 1.95 * eV, 1.97 * eV, 2.00 * eV, 2.12 * eV, 2.18 * eV, 2.23 * eV, 2.23 * eV, 2.29 * eV, 2.30 * eV, 2.37 * eV, 2.37 * eV, 2.39 * eV, 2.47 * eV, 2.47 * eV, 2.55 * eV, 2.62 * eV, 2.71 * eV, 2.80 * eV, 2.80 * eV, 2.94 * eV, 2.96 * eV, 2.99 * eV, 3.09 * eV, 3.17 * eV, 3.24 * eV, 3.30 * eV, 3.31 * eV, 3.36 * eV, 3.43 * eV, 3.44 * eV, 3.50 * eV, 3.53 * eV, 3.57 * eV, 3.58 * eV, 3.58 * eV, 3.64 * eV, 3.64 * eV, 3.65 * eV, 3.70 * eV, 3.72 * eV, 3.80 * eV, 3.82 * eV, 3.87 * eV, 3.91 * eV};
  G4double effi_mppc[62] = {0.03,    0.04,    0.05,    0.05,    0.06,    0.06,    0.07,    0.08,    0.08,    0.09,    0.11,    0.12,    0.12,    0.14,    0.15,    0.15,    0.19,    0.21,    0.22,    0.23,    0.24,    0.28,    0.30,    0.32,    0.32,    0.34,    0.34,    0.37,    0.36,    0.37,    0.39,    0.39,    0.39,    0.40,    0.40,    0.40,    0.39,    0.38,    0.38,    0.37,    0.35,    0.34,    0.31,    0.29,    0.30,    0.27,    0.25,    0.25,    0.22,    0.20,    0.21,    0.18,    0.18,    0.16,    0.14,    0.12,    0.12,    0.09,    0.05,    0.04,    0.02,    0.02};
  for (int i = 0; i < 62; i++)
    effi_mppc[i] *= _sipm_eff;

  G4MaterialPropertiesTable *table = fSiPMSurface3->GetMaterialPropertiesTable();
  if (table)
  {
    table->RemoveProperty("EFFICIENCY");
    table->AddProperty("EFFICIENCY", p_mppc, effi_mppc, 62);
  }
  else
  {
    table = new G4MaterialPropertiesTable();
    table->AddProperty("EFFICIENCY", p_mppc, effi_mppc, 62);
    fSiPMSurface3->SetMaterialPropertiesTable(table);
  }

  G4MaterialPropertiesTable *table4 = fSiPMSurface4->GetMaterialPropertiesTable();
  if (table4)
  {
    table4->RemoveProperty("EFFICIENCY");
    table4->AddProperty("EFFICIENCY", p_mppc, effi_mppc, 62);
  }
  else
  {
    table4 = new G4MaterialPropertiesTable();
    table4->AddProperty("EFFICIENCY", p_mppc, effi_mppc, 62);
    fSiPMSurface4->SetMaterialPropertiesTable(table4);
  }
}

void LYSimDetectorConstruction::SetTileAlpha(const double a)
{
  fTileBulkSurface->SetSigmaAlpha(a);
  _tile_alpha = a;
}

void LYSimDetectorConstruction::SetDimpleAlpha(const double a)
{
  fTileDimpleSurface->SetSigmaAlpha(a);
  _dimple_alpha = a;
}
