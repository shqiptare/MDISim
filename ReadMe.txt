MDISim/ReadMe.txt        see   make_MDISim_dir.sh

Tools for flexible optimisation of IR designs with application to FCC
Helmut Burkhardt / CERN, Manuela Boscolo / INFN-Frascati  
https://jacowfs.jlab.org/conf/proceedings/IPAC2015/papers/tupty031.pdf
MDISIm_dir is generated with make_MDISim_dir.sh

Works with standard root6     + extra MDI libs and headers loaded in root


Install on MAC
--------------
cp /afs/cern.ch/user/h/hbu/public/Manuela/MDISim_dir_mac.tgz   in your preferred location
tar -xzf MDISim_dir_mac.tgz
cd MDISim_dir ; source setup.sh


Install on lxplus
-----------------
cp /afs/cern.ch/user/h/hbu/public/Manuela/MDISim_dir_linux.tgz   in your preferred location
tar -xzf MDISim_dir_mac.tgz
cd MDISim_dir ; source setup.sh




Code files and their function
-----------------------------
CMakeLists.txt   steering CMake, also used to re-create MDISim.xcodeproj from time to time

main_, MyNtuple2Geom.C, .h   construct root geometry from MAD-X twiss and survey tfs
Twiss_util.C, Twiss_util.h   use twiss header and ntuple, beam parameters and SR

main_AddTracksToGeom.C       display tracks with geomettry in root
main_CLIC.C                  MDISim applied to CLIC
main_FCC.C                   MDISim applied to FCC
main_Fields_from_tfs.C       Extract fields from tfs files for easy input to Geant4
main_LEP.C                   MDISim applied to LEP
main_LHC.C                   MDISim applied to LHC
myEve.C, .h                  for plotting with TEve
myGraph3D.C, .h              root 3D track and cone plotting

Other utilities, 2d plots, to some extend historic 
--------------------------------------------------
Apertures.C,.h                 classwith aperture info
IPpos.C,.h                     index of IP, CalcAveBeam
main_, MyGeoManager.C,h        Used to try out different shapes
Plot.C, Plot.h                 2-d plots of elements
PlotAddTrack.C                 read/plot reads  TEveLine objects
makefile                       handmade makefile -  normally rather use CMake
main_My_VGM.C                  attempts to use VGM
main_, TwissTrack_Param.C, .h  2-d aperture and track plotting
myGraph.C,.h                   2-d aperture and track plotting
