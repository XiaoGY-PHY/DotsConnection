----------> uses
# use GaudiInterface GaudiInterface-01-* External
#   use GaudiKernel *  (no_version_directory)
#     use GaudiPolicy *  (no_version_directory)
#       use LCG_Settings *  (no_version_directory)
#       use Python * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#       use tcmalloc * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.7p3)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use libunwind v* LCG_Interfaces (no_version_directory) (native_version=5c2cade)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#     use Reflex * LCG_Interfaces (no_version_directory)
#       use LCG_Configuration v*  (no_version_directory)
#         use LCG_Platforms *  (no_version_directory)
#       use LCG_Settings v*  (no_version_directory)
#       use ROOT v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#     use Boost * LCG_Interfaces (no_version_directory) (native_version=1.50.0_python2.7)
#       use LCG_Configuration v*  (no_version_directory)
#       use LCG_Settings v*  (no_version_directory)
#       use Python v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#   use GaudiSvc *  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use CLHEP * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.0.4.5)
#     use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#     use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
# use BesCLHEP BesCLHEP-* External
#   use CLHEP v* LCG_Interfaces (no_version_directory) (native_version=2.0.4.5)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use HepPDT * LCG_Interfaces (no_version_directory) (native_version=2.06.01)
#     use LCG_Configuration v*  (no_version_directory)
#     use LCG_Settings v*  (no_version_directory)
#   use BesExternalArea BesExternalArea-* External
# use MdcGeomSvc MdcGeomSvc-* Mdc
#   use BesPolicy BesPolicy-01-* 
#     use BesCxxPolicy BesCxxPolicy-* 
#       use GaudiPolicy v*  (no_version_directory)
#     use BesFortranPolicy BesFortranPolicy-* 
#       use LCG_Settings v*  (no_version_directory)
#   use GaudiInterface GaudiInterface-* External
#   use calibUtil * Calibration
#     use GaudiInterface GaudiInterface-01-* External
#     use facilities * Calibration
#       use BesPolicy BesPolicy-* 
#     use xmlBase * Calibration
#       use BesPolicy * 
#       use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#       use facilities * Calibration
#     use rdbModel * Calibration
#       use BesPolicy * 
#       use facilities * Calibration
#       use xmlBase * Calibration
#       use MYSQL * External
#         use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#     use BesROOT BesROOT-00-* External
#       use CASTOR v* LCG_Interfaces (no_version_directory) (native_version=2.1.13-6)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#       use ROOT v* LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use GCCXML v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=0.9.0_20120309p2)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#         use Python v* LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#         use xrootd v* LCG_Interfaces (no_version_directory) (native_version=3.2.7)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#     use DatabaseSvc DatabaseSvc-* Database
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#       use sqlite * LCG_Interfaces (no_version_directory) (native_version=3070900)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#       use BesROOT * External
#   use CalibData * Calibration
#     use facilities facilities-* Calibration
#     use GaudiInterface * External
#     use BesROOT BesROOT-00-* External
#   use EventModel EventModel-* Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
# use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-01-* External
#   use CalibData CalibData-* Calibration
#   use CalibSvc CalibSvc-* Calibration
#     use CalibMySQLCnv * Calibration/CalibSvc
#       use BesPolicy * 
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use GaudiInterface * External
#       use MYSQL MYSQL-00-* External
#       use CalibDataSvc * Calibration/CalibSvc
#         use BesPolicy * 
#         use BesROOT * External
#         use calibUtil * Calibration
#         use CalibData * Calibration
#         use DstEvent DstEvent-* Event
#           use BesPolicy BesPolicy-* 
#           use GaudiInterface GaudiInterface-* External
#           use BesCLHEP BesCLHEP-* External
#           use EventModel EventModel-* Event
#         use EventModel EventModel-* Event
#         use GaudiKernel *  (no_version_directory)
#     use CalibROOTCnv * Calibration/CalibSvc
#       use BesPolicy * 
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use GaudiInterface * External
#       use BesROOT BesROOT-00-* External
#       use EventModel EventModel-* Event
#       use CalibDataSvc * Calibration/CalibSvc
#       use CalibMySQLCnv * Calibration/CalibSvc
#     use CalibDataSvc * Calibration/CalibSvc
#     use CalibTreeCnv * Calibration/CalibSvc
#       use BesPolicy * 
#       use calibUtil * Calibration
#       use CalibData * Calibration
#       use MYSQL MYSQL-00-* External
#       use GaudiInterface * External
#       use BesROOT BesROOT-00-* External
#       use DstEvent DstEvent-* Event
#       use EventModel EventModel-* Event
#       use CalibDataSvc * Calibration/CalibSvc
#       use CalibMySQLCnv * Calibration/CalibSvc
#       use DatabaseSvc DatabaseSvc-* Database
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use BesCLHEP BesCLHEP-* External
# use RawDataProviderSvc RawDataProviderSvc-* Event
#   use BesPolicy BesPolicy-01-* 
#   use RootPolicy RootPolicy-* 
#     use BesPolicy BesPolicy-* 
#     use BesROOT BesROOT-00-* External
#   use GaudiInterface GaudiInterface-01-* External
#   use BesCLHEP BesCLHEP-* External
#   use DetVerSvc DetVerSvc-* Utilities
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#   use MdcRawEvent MdcRawEvent-* Mdc
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawEvent RawEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use Identifier Identifier-* DetectorDescription
#         use BesPolicy BesPolicy-* 
#       use EventModel EventModel-* Event
#     use EventModel EventModel-* Event
#   use MdcRecEvent MdcRecEvent-* Mdc
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcGeomSvc MdcGeomSvc-* Mdc
#     use RelTable RelTable-* Event
#       use BesPolicy BesPolicy-01-* 
#       use GaudiInterface GaudiInterface-01-* External
#     use EventModel EventModel-* Event
#     use Identifier Identifier-* DetectorDescription
#     use DstEvent DstEvent-* Event
#   use TofRawEvent TofRawEvent-* Tof
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawEvent RawEvent-* Event
#     use EventModel EventModel-* Event
#   use TofCaliSvc TofCaliSvc-* Tof
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use CalibSvc * Calibration
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use MYSQL MYSQL-00-00-* External
#     use BesROOT BesROOT-* External
#   use TofQCorrSvc TofQCorrSvc-* Tof
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use CalibData * Calibration
#     use EventModel EventModel-* Event
#     use DatabaseSvc DatabaseSvc-* Database
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use MYSQL MYSQL-00-00-* External
#     use BesROOT BesROOT-* External
#   use TofQElecSvc TofQElecSvc-* Tof
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use calibUtil * Calibration
#     use CalibData * Calibration
#     use CalibSvc * Calibration
#     use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#     use MYSQL MYSQL-00-00-* External
#     use BesROOT BesROOT-* External
#   use EmcRawEvent EmcRawEvent-* Emc
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawEvent RawEvent-* Event
#     use EventModel EventModel-* Event
#     use Identifier Identifier-* DetectorDescription
#     use EmcWaveform EmcWaveform-* Emc
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#   use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-* External
#     use CalibData CalibData-* Calibration
#     use CalibDataSvc CalibDataSvc-* Calibration/CalibSvc
#     use CalibROOTCnv CalibROOTCnv-* Calibration/CalibSvc
#     use EmcRecGeoSvc EmcRecGeoSvc-* Emc
#       use BesPolicy BesPolicy-* 
#       use Identifier Identifier-* DetectorDescription
#       use ROOTGeo ROOTGeo-* DetectorDescription
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use BesCLHEP BesCLHEP-* External
#         use BesROOT BesROOT-* External
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use GdmlToRoot GdmlToRoot-* External
#           use BesExternalArea BesExternalArea-* External
#           use BesROOT BesROOT-* External
#         use GdmlManagement GdmlManagement-* DetectorDescription
#           use BesExternalArea BesExternalArea-* External
#       use BesCLHEP BesCLHEP-* External
#       use GaudiInterface GaudiInterface-* External
#     use EmcGeneralClass EmcGeneralClass-* Emc
#       use BesPolicy BesPolicy-* 
#       use Identifier Identifier-* DetectorDescription
#     use BesCLHEP BesCLHEP-* External
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use EvTimeEvent EvTimeEvent-* Event
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use BesCLHEP BesCLHEP-* External
#     use MdcGeomSvc MdcGeomSvc-* Mdc
#     use RelTable RelTable-* Event
#     use EventModel EventModel-* Event
#     use Identifier Identifier-* DetectorDescription
# use MdcUtilitySvc MdcUtilitySvc-* Mdc/MdcCheckUtil
#   use OfflineEventLoopMgr * Control
#     use BesPolicy BesPolicy-* 
#     use BesBoost BesBoost-* External
#       use Boost v* LCG_Interfaces (no_version_directory) (native_version=1.50.0_python2.7)
#     use GaudiInterface GaudiInterface-01-* External
#     use GaudiCoreSvc v*  (no_version_directory)
#       use GaudiKernel *  (no_version_directory)
#       use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#     use GaudiCommonSvc v*  (no_version_directory)
#       use GaudiKernel *  (no_version_directory)
#       use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#       use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#       use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#     use EventModel EventModel-* Event
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-* External
#   use BesCLHEP BesCLHEP-* External
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use RawDataCnv RawDataCnv-* Event
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use RawFile RawFile-* Event
#       use BesPolicy BesPolicy-* 
#       use BesExternalArea BesExternalArea-* External
#       use IRawFile IRawFile-* Event
#         use BesPolicy BesPolicy-01-* 
#       use XmlRpcC XmlRpcC-* External
#         use BesExternalArea BesExternalArea-* External
#       use mysql * LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#     use EventModel EventModel-* Event
#     use RawDataCnvBase RawDataCnvBase-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use EventModel EventModel-* Event
#     use RawEvent RawEvent-00-* Event
#     use McTruth McTruth-* Event
#       use BesPolicy BesPolicy-01-* 
#       use EventModel EventModel-* Event
#       use GaudiInterface GaudiInterface-01-* External
#       use Identifier Identifier-* DetectorDescription
#       use RelTable RelTable-* Event
#     use Identifier Identifier-* DetectorDescription
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use MucRawEvent MucRawEvent-* Muc
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use RawEvent RawEvent-* Event
#       use EventModel EventModel-* Event
#     use TrigEvent TrigEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#     use HltEvent HltEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use RawEvent RawEvent-* Event
#       use EventModel EventModel-* Event
#       use Identifier Identifier-* DetectorDescription
#     use LumiDigi LumiDigi-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-01-* External
#       use RawEvent RawEvent-* Event
#       use TofRawEvent TofRawEvent-* Tof
#     use ZddEvent ZddEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#     use eformat eformat-* Event
#       use BesPolicy BesPolicy-* 
#       use ers ers-* Event
#         use BesPolicy BesPolicy-* 
#     use TofSim TofSim-* Simulation/BOOST
#       use BesPolicy BesPolicy-01-* 
#       use BesGeant4 BesGeant4-00-* External
#         use BesExternalArea BesExternalArea-00-* External
#         use BesCLHEP BesCLHEP-00-* External
#       use GdmlToG4 GdmlToG4-* External
#         use BesExternalArea BesExternalArea-* External
#         use BesGeant4 BesGeant4-* External
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#       use SimUtil SimUtil-* Simulation/BOOST
#         use BesPolicy BesPolicy-01-* 
#         use BesGeant4 BesGeant4-00-* External
#       use TruSim TruSim-* Simulation/BOOST
#         use BesPolicy BesPolicy-01-* 
#         use BesGeant4 BesGeant4-* External
#         use BesCLHEP BesCLHEP-* External
#         use GeneratorObject GeneratorObject-* Generator
#           use BesPolicy BesPolicy-* 
#           use HepMC * LCG_Interfaces (no_version_directory) (native_version=2.06.08)
#           use EventModel EventModel-* Event
#       use GaudiInterface GaudiInterface-* External
#       use G4Svc G4Svc-* Simulation
#         use BesCLHEP BesCLHEP-* External
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use GeneratorObject GeneratorObject-* Generator
#         use BesGeant4 BesGeant4-00-* External
#         use RealizationSvc RealizationSvc-* Simulation/Realization
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-01-* External
#           use BesCLHEP * External
#           use EmcCalibConstSvc EmcCalibConstSvc-* Emc
#           use EventModel EventModel-* Event
#           use DetVerSvc DetVerSvc-* Utilities
#       use GdmlManagement GdmlManagement-* DetectorDescription
#       use G4Geo G4Geo-* DetectorDescription
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use BesCLHEP BesCLHEP-* External
#         use BesGeant4 BesGeant4-* External
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use GdmlToG4 GdmlToG4-* External
#         use GdmlManagement GdmlManagement-* DetectorDescription
#         use Identifier Identifier-* DetectorDescription
#         use SimUtil SimUtil-* Simulation/BOOST
#       use TofCaliSvc TofCaliSvc-* Tof
#       use TofSimSvc TofSimSvc-* Tof
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-* External
#         use calibUtil * Calibration
#         use CalibData * Calibration
#         use CalibSvc * Calibration
#         use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#         use MYSQL MYSQL-00-00-* External
#       use TofQElecSvc TofQElecSvc-* Tof
#     use RawDataProviderSvc RawDataProviderSvc-* Event
#   use McTruth McTruth-* Event
#   use RawDataProviderSvc RawDataProviderSvc-* Event
#   use RootIO RootIO-* Event
#     use GaudiKernel *  (no_version_directory)
#     use Gaudi *  (no_version_directory)
#       use GaudiSys *  (no_version_directory)
#         use GaudiKernel *  (no_version_directory)
#         use GaudiCoreSvc *  (no_auto_imports) (no_version_directory)
#         use GaudiCommonSvc *  (no_auto_imports) (no_version_directory)
#         use GaudiAud *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#         use GaudiAlg *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#           use GaudiUtils *  (no_version_directory)
#             use GaudiKernel *  (no_version_directory)
#             use ROOT * LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#             use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#             use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#             use uuid * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.42)
#               use LCG_Configuration v*  (no_version_directory)
#               use LCG_Settings v*  (no_version_directory)
#             use Reflex * LCG_Interfaces (no_auto_imports) (no_version_directory)
#             use XercesC * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.1.1p1)
#           use ROOT * LCG_Interfaces (no_version_directory) (native_version=5.34.09)
#           use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#           use Boost * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.50.0_python2.7)
#         use GaudiPython *  (no_auto_imports) (no_version_directory)
#           use GaudiKernel *  (no_version_directory)
#           use GaudiAlg *  (no_version_directory)
#           use GaudiUtils *  (no_version_directory)
#           use AIDA * LCG_Interfaces (no_version_directory) (native_version=3.2.1)
#           use CLHEP * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.0.4.5)
#           use Python * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=2.7.3)
#           use Reflex * LCG_Interfaces (no_auto_imports) (no_version_directory)
#       use pytools * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=1.8_python2.7)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use Python v* LCG_Interfaces (no_version_directory) (native_version=2.7.3)
#         use mysql v* LCG_Interfaces (no_version_directory) (native_version=5.5.14)
#     use ReconEvent ReconEvent-* Event
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#       use BesCLHEP BesCLHEP-* External
#       use EventModel EventModel-* Event
#       use ExtEvent ExtEvent-* Event
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use BesCLHEP BesCLHEP-* External
#         use EventModel EventModel-* Event
#         use DstEvent DstEvent-* Event
#     use RootCnvSvc RootCnvSvc-* Event
#       use BesCLHEP BesCLHEP-* External
#       use BesPolicy BesPolicy-01-* 
#       use DataInfoSvc DataInfoSvc-* Control
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#       use TagFilterSvc TagFilterSvc-* Event
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use BesROOT BesROOT-00-* External
#       use BesROOT BesROOT-00-* External
#       use GaudiInterface GaudiInterface-* External
#       use Identifier Identifier-* DetectorDescription
#       use RawEvent RawEvent-* Event
#       use RootEventData RootEventData-* Event
#         use RootPolicy RootPolicy-* 
#         use BesROOT BesROOT-* External
#         use MucRecEvent MucRecEvent-* Muc
#           use BesPolicy BesPolicy-01-* 
#           use GaudiInterface GaudiInterface-01-* External
#           use Identifier Identifier-* DetectorDescription
#           use EventModel EventModel-* Event
#           use ExtEvent ExtEvent-* Event
#           use MucGeomSvc MucGeomSvc-* Muc
#             use BesPolicy BesPolicy-01-* 
#             use GaudiInterface GaudiInterface-* External
#             use Identifier Identifier-* DetectorDescription
#             use ROOTGeo ROOTGeo-* DetectorDescription
#             use BesCLHEP BesCLHEP-* External
#             use BesROOT BesROOT-* External
#             use XercesC * LCG_Interfaces (no_version_directory) (native_version=3.1.1p1)
#             use GdmlToRoot GdmlToRoot-* External
#             use G4Geo G4Geo-* DetectorDescription
#           use DstEvent DstEvent-* Event
#         use Identifier Identifier-* DetectorDescription
#       use MdcRawEvent MdcRawEvent-* Mdc
#       use EventModel EventModel-* Event
#       use EmcRawEvent EmcRawEvent-* Emc
#       use TofRawEvent TofRawEvent-* Tof
#       use TofRecEvent TofRecEvent-* Tof
#         use BesPolicy BesPolicy-01-* 
#         use GaudiInterface GaudiInterface-01-* External
#         use Identifier Identifier-* DetectorDescription
#         use EventModel EventModel-* Event
#         use DstEvent * Event
#       use MucRawEvent MucRawEvent-* Muc
#       use MucRecEvent MucRecEvent-* Muc
#       use EmcRecEventModel EmcRecEventModel-* Emc
#         use BesPolicy BesPolicy-* 
#         use Identifier Identifier-* DetectorDescription
#         use BesCLHEP BesCLHEP-* External
#         use EventModel EventModel-* Event
#         use DstEvent DstEvent-* Event
#         use EmcRecGeoSvc EmcRecGeoSvc-* Emc
#       use MdcRecEvent MdcRecEvent-* Mdc
#       use EvTimeEvent EvTimeEvent-* Event
#       use ZddEvent ZddEvent-* Event
#       use EvtRecEvent EvtRecEvent-* Event
#         use BesPolicy BesPolicy-* 
#         use GaudiInterface GaudiInterface-* External
#         use BesCLHEP BesCLHEP-* External
#         use EventModel EventModel-* Event
#         use EvTimeEvent EvTimeEvent-* Event
#         use MdcRecEvent MdcRecEvent-* Mdc
#         use TofRecEvent TofRecEvent-* Tof
#         use EmcRecEventModel EmcRecEventModel-* Emc
#         use MucRecEvent MucRecEvent-* Muc
#         use ExtEvent ExtEvent-* Event
#         use DstEvent DstEvent-* Event
#       use EventNavigator EventNavigator-* Event
#         use BesPolicy * 
#         use GaudiInterface * External
#         use McTruth McTruth-* Event
#         use EmcRecEventModel * Emc
#         use MdcRecEvent * Mdc
#         use MdcRawEvent * Mdc
#         use MucRecEvent * Muc
#         use TofRecEvent * Tof
#       use McTruth McTruth-* Event
#       use ExtEvent ExtEvent-* Event
#       use DstEvent DstEvent-* Event
#       use ReconEvent ReconEvent-* Event
#       use TrigEvent TrigEvent-* Event
#       use HltEvent HltEvent-* Event
#       use LumiDigi LumiDigi-* Event
#     use BesROOT BesROOT-00-* External
#     use BesCLHEP BesCLHEP-* External
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use TofRawEvent TofRawEvent-* Tof
#     use EmcRawEvent EmcRawEvent-* Emc
#     use DstEvent DstEvent-* Event
#     use EventModel EventModel-* Event
#     use McTruth McTruth-* Event
#     use Identifier Identifier-* DetectorDescription
#     use RootPolicy RootPolicy-* 
#     use RootEventData RootEventData-* Event
#   use MagneticField MagneticField-* 
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use BesCLHEP * External
#     use BesROOT * External
#     use EventModel EventModel-* Event
#     use rdbModel * Calibration
#     use DatabaseSvc * Database
#     use BesTimerSvc BesTimerSvc-* Utilities
#       use BesPolicy BesPolicy-* 
#       use GaudiInterface GaudiInterface-* External
#   use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use Identifier Identifier-* DetectorDescription
#     use MdcGeomSvc MdcGeomSvc-* Mdc
#     use BesCLHEP BesCLHEP-* External
#   use TrkBase TrkBase-* Reconstruction/MdcPatRec
#     use BesPolicy BesPolicy-01-* 
#     use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#     use MdcRecoUtil MdcRecoUtil-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-* 
#       use BesCLHEP BesCLHEP-* External
#     use ProxyDict ProxyDict-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-01-* 
#     use ProbTools ProbTools-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-01-* 
#     use BField BField-* Reconstruction/MdcPatRec
#       use BesPolicy BesPolicy-* 
#       use MdcGeom MdcGeom-* Reconstruction/MdcPatRec
#       use MagneticField MagneticField-* 
#     use BesCLHEP BesCLHEP-* External
#     use CERNLIB CERNLIB-* External
#       use cernlib v* LCG_Interfaces (no_version_directory) (native_version=2006a)
#         use LCG_Configuration v*  (no_version_directory)
#         use LCG_Settings v*  (no_version_directory)
#         use blas v* LCG_Interfaces (no_version_directory) (native_version=20110419)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#         use lapack v* LCG_Interfaces (no_version_directory) (native_version=3.4.0)
#           use LCG_Configuration v*  (no_version_directory)
#           use LCG_Settings v*  (no_version_directory)
#           use blas v* LCG_Interfaces (no_version_directory) (native_version=20110419)
#       use CASTOR v* LCG_Interfaces (no_version_directory) (native_version=2.1.13-6)
#       use Xt * External (native_version=X11R6)
#     use BesBoost BesBoost-* External
#     use MdcRecEvent MdcRecEvent-* Mdc
#   use TrackUtil TrackUtil-* Reconstruction
#     use MagneticField MagneticField-* 
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use BesCLHEP BesCLHEP-* External
#     use MdcRawEvent MdcRawEvent-* Mdc
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use CERNLIB CERNLIB-* External
# use MdcRawEvent MdcRawEvent-* Mdc
# use KalFitAlg KalFitAlg-* Reconstruction
#   use BesPolicy BesPolicy-01-* 
#   use GaudiInterface GaudiInterface-* External
#   use CERNLIB CERNLIB-* External
#   use MdcRawEvent MdcRawEvent-* Mdc
#   use BesCLHEP BesCLHEP-* External
#   use EventModel EventModel-* Event
#   use EvTimeEvent EvTimeEvent-* Event
#   use MdcTables MdcTables-* Mdc
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use BesPolicy BesPolicy-* 
#     use GaudiInterface GaudiInterface-01-* External
#     use MdcGeomSvc MdcGeomSvc-* Mdc
#     use Identifier Identifier-* DetectorDescription
#   use ReconEvent ReconEvent-* Event
#   use McTruth McTruth-* Event
#   use G4Svc G4Svc-00-* Simulation
#   use G4Geo * DetectorDescription
#   use RootHistCnv v*  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use AIDA * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=3.2.1)
#     use ROOT * LCG_Interfaces (no_auto_imports) (no_version_directory) (native_version=5.34.09)
#   use BesROOT BesROOT-* External
#   use BesAIDA * External
#     use AIDA v* LCG_Interfaces (no_version_directory) (native_version=3.2.1)
#       use LCG_Configuration v*  (no_version_directory)
#       use LCG_Settings v*  (no_version_directory)
#   use MdcRecEvent MdcRecEvent-* Mdc
#   use MdcGeomSvc MdcGeomSvc-* Mdc
#   use MdcCalibFunSvc MdcCalibFunSvc-* Mdc
#   use MagneticField MagneticField-* 
#   use PartPropSvc *  (no_version_directory)
#     use GaudiPolicy *  (no_version_directory)
#     use GaudiKernel *  (no_version_directory)
#     use HepPDT * LCG_Interfaces (no_version_directory) (native_version=2.06.01)
#   use VertexFit VertexFit-* Analysis
#     use BesPolicy BesPolicy-01-* 
#     use GaudiInterface GaudiInterface-* External
#     use BesCLHEP BesCLHEP-* External
#     use MYSQL MYSQL-* External
#     use DstEvent DstEvent-* Event
#     use MdcRecEvent MdcRecEvent-* Mdc
#     use EmcRecEventModel EmcRecEventModel-* Emc
#     use MagneticField MagneticField-* 
#     use EventModel EventModel-* Event
#     use EvtRecEvent EvtRecEvent-* Event
#     use DatabaseSvc DatabaseSvc-* Database
#
# Selection :
use CMT v1r25 (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib)
use Xt Xt-00-00-02 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use BesExternalArea BesExternalArea-00-00-22 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use XmlRpcC XmlRpcC-00-00-02 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use GdmlManagement GdmlManagement-00-00-44 DetectorDescription (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use LCG_Platforms v1  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use LCG_Configuration v1  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use LCG_Settings v1  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use blas v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use lapack v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use cernlib v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use uuid v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use AIDA v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesAIDA BesAIDA-00-00-01 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use sqlite v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use CASTOR v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use CERNLIB CERNLIB-01-02-03 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use mysql v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use MYSQL MYSQL-00-00-09 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use XercesC v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesFortranPolicy BesFortranPolicy-00-01-03  (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use HepPDT v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use HepMC v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use CLHEP v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesCLHEP BesCLHEP-00-00-11 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use BesGeant4 BesGeant4-00-00-11 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use GdmlToG4 GdmlToG4-00-00-11 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use xrootd v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use GCCXML v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use libunwind v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use tcmalloc v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use Python v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use pytools v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a) (no_auto_imports)
use Boost v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesBoost BesBoost-00-00-01 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ROOT v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use BesROOT BesROOT-00-00-08 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use GdmlToRoot GdmlToRoot-00-00-14 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use Reflex v1 LCG_Interfaces (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a)
use GaudiPolicy v12r7  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use BesCxxPolicy BesCxxPolicy-00-01-01  (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use BesPolicy BesPolicy-01-05-05  (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ProbTools ProbTools-00-00-01 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ProxyDict ProxyDict-00-00-01 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcRecoUtil MdcRecoUtil-00-01-08 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use SimUtil SimUtil-00-00-37 Simulation/BOOST (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ers ers-00-00-03 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use eformat eformat-00-00-04 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use IRawFile IRawFile-00-00-05 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RawFile RawFile-00-00-14 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use Identifier Identifier-00-02-17 DetectorDescription (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcGeneralClass EmcGeneralClass-00-00-04 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RootPolicy RootPolicy-00-01-03  (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use facilities facilities-00-00-04 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use xmlBase xmlBase-00-00-03 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use rdbModel rdbModel-00-01-01 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use GaudiKernel v28r8  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use PartPropSvc v4r6  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use RootHistCnv v11r2  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiUtils v4r6  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiAlg v14r6  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiPython v12r6  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiAud v9r9  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9) (no_auto_imports)
use GaudiCommonSvc v1r6  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiCoreSvc v1r5  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiSys v23r9  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use Gaudi v23r9  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiSvc v19r4  (/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9)
use GaudiInterface GaudiInterface-01-03-07 External (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use BesTimerSvc BesTimerSvc-00-00-12 Utilities (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TagFilterSvc TagFilterSvc-00-00-21 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use DataInfoSvc DataInfoSvc-00-00-03 Control (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use G4Geo G4Geo-00-00-13 DetectorDescription (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ZddEvent ZddEvent-00-00-04 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TrigEvent TrigEvent-00-01-02 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ROOTGeo ROOTGeo-00-00-15 DetectorDescription (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MucGeomSvc MucGeomSvc-00-02-25 Muc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcRecGeoSvc EmcRecGeoSvc-01-01-07 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcWaveform EmcWaveform-00-00-03 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RelTable RelTable-00-00-02 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use DetVerSvc DetVerSvc-00-00-07 Utilities (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EventModel EventModel-01-05-36 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use GeneratorObject GeneratorObject-00-01-05 Generator (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TruSim TruSim-00-00-17 Simulation/BOOST (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use McTruth McTruth-00-02-21 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RawDataCnvBase RawDataCnvBase-01-00-03 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use OfflineEventLoopMgr OfflineEventLoopMgr-00-00-15 Control (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RawEvent RawEvent-00-03-19 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use HltEvent HltEvent-00-02-07 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MucRawEvent MucRawEvent-00-02-02 Muc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcRawEvent EmcRawEvent-00-02-06 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofRawEvent TofRawEvent-00-02-07 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use LumiDigi LumiDigi-00-00-02 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcRawEvent MdcRawEvent-00-03-08 Mdc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use DstEvent DstEvent-00-02-51 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcRecEventModel EmcRecEventModel-01-01-18 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofRecEvent TofRecEvent-00-02-14 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ExtEvent ExtEvent-00-00-13 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MucRecEvent MucRecEvent-00-02-52 Muc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RootEventData RootEventData-00-03-85 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use ReconEvent ReconEvent-00-00-04 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibData CalibData-00-01-24 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use DatabaseSvc DatabaseSvc-00-00-26 Database (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MagneticField MagneticField-00-02-04  (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofQCorrSvc TofQCorrSvc-00-00-10 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use calibUtil calibUtil-00-00-46 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibDataSvc CalibDataSvc-00-01-07 Calibration/CalibSvc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibMySQLCnv CalibMySQLCnv-00-02-05 Calibration/CalibSvc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibTreeCnv CalibTreeCnv-00-02-02 Calibration/CalibSvc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibROOTCnv CalibROOTCnv-00-02-03 Calibration/CalibSvc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EmcCalibConstSvc EmcCalibConstSvc-00-00-13 Emc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RealizationSvc RealizationSvc-00-00-48 Simulation/Realization (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use G4Svc G4Svc-00-01-53 Simulation (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use CalibSvc CalibSvc-00-03-04 Calibration (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofSimSvc TofSimSvc-00-00-04 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofQElecSvc TofQElecSvc-00-00-05 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofCaliSvc TofCaliSvc-00-01-19 Tof (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TofSim TofSim-00-02-37 Simulation/BOOST (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcGeomSvc MdcGeomSvc-00-01-37 Mdc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcGeom MdcGeom-00-01-17 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use BField BField-00-01-02 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EvTimeEvent EvTimeEvent-00-00-08 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcRecEvent MdcRecEvent-00-05-14 Mdc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcTables MdcTables-00-00-11 Mdc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TrackUtil TrackUtil-00-00-08 Reconstruction (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use TrkBase TrkBase-00-01-12 Reconstruction/MdcPatRec (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EventNavigator EventNavigator-00-01-03 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use EvtRecEvent EvtRecEvent-00-02-03 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use VertexFit VertexFit-00-02-87 Analysis (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RootCnvSvc RootCnvSvc-03-00-06 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RootIO RootIO-00-01-34 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcCalibFunSvc MdcCalibFunSvc-00-03-16 Mdc (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use KalFitAlg KalFitAlg-00-07-55-p03 Reconstruction (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RawDataProviderSvc RawDataProviderSvc-00-03-50 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use RawDataCnv RawDataCnv-00-05-11 Event (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
use MdcUtilitySvc MdcUtilitySvc-00-00-09 Mdc/MdcCheckUtil (/cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a)
----------> tags
CMTv1 (from CMTVERSION)
CMTr25 (from CMTVERSION)
CMTp0 (from CMTVERSION)
Linux (from uname) package [CMT LCG_Platforms BesPolicy] implies [Unix host-linux]
x86_64-slc6-gcc46-opt (from CMTCONFIG) package [LCG_Platforms] implies [target-x86_64 target-slc6 target-gcc46 target-opt]
mg20220135_no_config (from PROJECT) excludes [mg20220135_config]
mg20220135_root (from PROJECT) excludes [mg20220135_no_root]
mg20220135_cleanup (from PROJECT) excludes [mg20220135_no_cleanup]
mg20220135_scripts (from PROJECT) excludes [mg20220135_no_scripts]
mg20220135_prototypes (from PROJECT) excludes [mg20220135_no_prototypes]
mg20220135_with_installarea (from PROJECT) excludes [mg20220135_without_installarea]
mg20220135_with_version_directory (from PROJECT) excludes [mg20220135_without_version_directory]
mg20220135 (from PROJECT)
BOSS_no_config (from PROJECT) excludes [BOSS_config]
BOSS_root (from PROJECT) excludes [BOSS_no_root]
BOSS_cleanup (from PROJECT) excludes [BOSS_no_cleanup]
BOSS_scripts (from PROJECT) excludes [BOSS_no_scripts]
BOSS_no_prototypes (from PROJECT) excludes [BOSS_prototypes]
BOSS_with_installarea (from PROJECT) excludes [BOSS_without_installarea]
BOSS_with_version_directory (from PROJECT) excludes [BOSS_without_version_directory]
GAUDI_no_config (from PROJECT) excludes [GAUDI_config]
GAUDI_root (from PROJECT) excludes [GAUDI_no_root]
GAUDI_cleanup (from PROJECT) excludes [GAUDI_no_cleanup]
GAUDI_scripts (from PROJECT) excludes [GAUDI_no_scripts]
GAUDI_prototypes (from PROJECT) excludes [GAUDI_no_prototypes]
GAUDI_with_installarea (from PROJECT) excludes [GAUDI_without_installarea]
GAUDI_without_version_directory (from PROJECT) excludes [GAUDI_with_version_directory]
LCGCMT_no_config (from PROJECT) excludes [LCGCMT_config]
LCGCMT_no_root (from PROJECT) excludes [LCGCMT_root]
LCGCMT_cleanup (from PROJECT) excludes [LCGCMT_no_cleanup]
LCGCMT_scripts (from PROJECT) excludes [LCGCMT_no_scripts]
LCGCMT_prototypes (from PROJECT) excludes [LCGCMT_no_prototypes]
LCGCMT_without_installarea (from PROJECT) excludes [LCGCMT_with_installarea]
LCGCMT_with_version_directory (from PROJECT) excludes [LCGCMT_without_version_directory]
x86_64 (from package CMT) package [LCG_Platforms] implies [host-x86_64]
slc78 (from package CMT)
gcc463 (from package CMT)
Unix (from package CMT) package [LCG_Platforms] implies [host-unix] excludes [WIN32 Win32]
c_native_dependencies (from package CMT) activated GaudiPolicy
cpp_native_dependencies (from package CMT) activated GaudiPolicy
target-unix (from package LCG_Settings) activated LCG_Platforms
target-x86_64 (from package LCG_Settings) activated LCG_Platforms
target-gcc46 (from package LCG_Settings) package [LCG_Platforms] implies [target-gcc4 target-lcg-compiler lcg-compiler] activated LCG_Platforms
host-x86_64 (from package LCG_Settings) activated LCG_Platforms
target-slc (from package LCG_Settings) package [LCG_Platforms] implies [target-linux] activated LCG_Platforms
target-gcc (from package LCG_Settings) activated LCG_Platforms
target-gcc4 (from package LCG_Settings) package [LCG_Platforms] implies [target-gcc] activated LCG_Platforms
target-lcg-compiler (from package LCG_Settings) activated LCG_Platforms
host-linux (from package LCG_Platforms) package [LCG_Platforms] implies [host-unix]
host-unix (from package LCG_Platforms)
target-opt (from package LCG_Platforms)
target-slc6 (from package LCG_Platforms) package [LCG_Platforms] implies [target-slc]
target-linux (from package LCG_Platforms) package [LCG_Platforms] implies [target-unix]
lcg-compiler (from package LCG_Platforms)
ROOT_GE_5_15 (from package LCG_Configuration)
ROOT_GE_5_19 (from package LCG_Configuration)
HasAthenaRunTime (from package BesPolicy)
ROOTBasicLibs (from package BesROOT)
----------> CMTPATH
# Add path /workfs2/bes/mg20220135/Boss706a from initialization
# Add path /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.6.a from initialization
# Add path /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/gaudi/GAUDI_v23r9 from initialization
# Add path /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/ExternalLib/LCGCMT/LCGCMT_65a from initialization
