# Microsoft Developer Studio Project File - Name="Diffusion" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Diffusion - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Diffusion.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Diffusion.mak" CFG="Diffusion - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Diffusion - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "Diffusion - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Diffusion - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /include:"Release/" /compile_only /nologo /warn:nofileopt
# ADD F90 /real_size:64 /define:"VF" /include:"Release/" /compile_only /nologo /check:bounds /warn:nofileopt
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /stack:0x1e8480 /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /warn:nofileopt
# ADD F90 /real_size:64 /browser /define:"VF" /include:"Debug/" /compile_only /nologo /debug:full /optimize:0 /check:bounds /warn:nofileopt
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Diffusion - Win32 Release"
# Name "Diffusion - Win32 Debug"
# Begin Group "Libs"

# PROP Default_Filter ""
# Begin Source File

SOURCE="C:\Program Files\DevStudio\DF\IMSL\LIB\SF90MP.LIB"
# End Source File
# Begin Source File

SOURCE="C:\Program Files\DevStudio\DF\IMSL\LIB\SMATHD.LIB"
# End Source File
# End Group
# Begin Source File

SOURCE=.\DoInX.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

DEP_F90_DOINX=\
	".\Release\isDiagPredomFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_DOINX=\
	".\Release\mpif.h"\
	

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

DEP_F90_DOINX=\
	".\Debug\isDiagPredomFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_DOINX=\
	".\Debug\mpif.h"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\DoInY.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

DEP_F90_DOINY=\
	".\Release\isDiagPredomFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_DOINY=\
	".\Release\mpif.h"\
	

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

DEP_F90_DOINY=\
	".\Debug\isDiagPredomFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_DOINY=\
	".\Debug\mpif.h"\
	

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\fRHS.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\GetMass.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\isDiagPredom.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\kron.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\SaveToDisk.f90
# End Source File
# Begin Source File

SOURCE=.\SwapXY.f90
# End Source File
# Begin Source File

SOURCE=.\SwapYX.f90
# End Source File
# Begin Source File

SOURCE=.\Test.f90

!IF  "$(CFG)" == "Diffusion - Win32 Release"

DEP_F90_TEST_=\
	".\Release\fRHSFunction.mod"\
	".\Release\GetMassSubroutine.mod"\
	".\Release\kronFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_TEST_=\
	".\Release\mpif.h"\
	

!ELSEIF  "$(CFG)" == "Diffusion - Win32 Debug"

DEP_F90_TEST_=\
	".\Debug\fRHSFunction.mod"\
	".\Debug\GetMassSubroutine.mod"\
	".\Debug\kronFunction.mod"\
	"C:\Program Files\DevStudio\DF\IMSL\Include\linear_operators.mod"\
	
NODEP_F90_TEST_=\
	".\Debug\mpif.h"\
	

!ENDIF 

# End Source File
# End Target
# End Project
