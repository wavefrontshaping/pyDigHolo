import ctypes as ct
import numpy as np
import matplotlib.pyplot as plt
import os
from typing import List, Union

DIGHOLO_ERROR_SUCCESS = 0 # Success (no error)
DIGHOLO_ERROR_ERROR = 1 #!< Generic error 
DIGHOLO_ERROR_INVALIDHANDLE = 2 #!< handleIdx provided doesn't exist 
DIGHOLO_ERROR_NULLPOINTER = 3 #!< A pointer provided was null 
DIGHOLO_ERROR_SETFRAMEBUFFERDISABLED = 4 #!< Requested to update the frame buffer, but it is currently disabled for update. Deprecated code.
DIGHOLO_ERROR_INVALIDDIMENSION = 5 #!< Dimension provided is invalid. (e.g. less than zero or greater than max possible value) 
DIGHOLO_ERROR_INVALIDPOLARISATION = 6 #!< Specified polarisation component doesn't exist. (e.g. there is only 1 polarisation component, but attempt to address second component) 
DIGHOLO_ERROR_INVALIDAXIS = 7 #!< Specified axis doesn't exist. (e.g. attempt to address third axis of a 2D system, or a negative axis) 
DIGHOLO_ERROR_INVALIDARGUMENT = 8  #!< A specified argument is not valid. (e.g. attempt to set value to out of range or otherwise meaningless value) 
DIGHOLO_ERROR_MEMORYALLOCATION = 9 #!< Memory allocation failed 
DIGHOLO_ERROR_FILENOTCREATED = 10 #!< File could not be created. (e.g. trying to redirect the console to a file path that doesn't exist)
DIGHOLO_ERROR_FILENOTFOUND = 11 #!< File could not be found/opened. (e.g. specifying a settings file path to open that doesn't exist)


class DigHoloError(Exception):
    pass
    
def check_error(ret):
    if ret == DIGHOLO_ERROR_ERROR:
        raise DigHoloError('Unknown error')
    elif ret == DIGHOLO_ERROR_INVALIDHANDLE:
        raise DigHoloError("HandleIdx provided invalid")
    elif ret == DIGHOLO_ERROR_NULLPOINTER:
        raise DigHoloError("A pointer provided was null")
    elif ret == DIGHOLO_ERROR_SETFRAMEBUFFERDISABLED:
        raise DigHoloError("Requested to update the frame buffer, but it is currently disabled for update. Deprecated code.")
    elif ret == DIGHOLO_ERROR_INVALIDDIMENSION:
        raise DigHoloError("Dimension provided is invalid. (e.g. less than zero or greater than max possible value) ")
    elif ret == DIGHOLO_ERROR_INVALIDPOLARISATION:
        raise DigHoloError("Specified polarisation component doesn't exist. (e.g. there is only 1 polarisation component, but attempt to address second component)")
    elif ret == DIGHOLO_ERROR_INVALIDAXIS:
        raise DigHoloError("Specified axis doesn't exist. (e.g. attempt to address third axis of a 2D system, or a negative axis) ")
    elif ret == DIGHOLO_ERROR_INVALIDARGUMENT:
        raise DigHoloError("A specified argument is not valid. (e.g. attempt to set value to out of range or otherwise meaningless value) ")
    elif ret == DIGHOLO_ERROR_MEMORYALLOCATION:
        raise DigHoloError("Memory allocation failed ")
    elif ret == DIGHOLO_ERROR_FILENOTCREATED:
        raise DigHoloError("File could not be created. (e.g. trying to redirect the console to a file path that doesn't exist)d")
    elif ret == DIGHOLO_ERROR_FILENOTFOUND:
        raise DigHoloError("File could not be found/opened. (e.g. specifying a settings file path to open that doesn't exist)")

class digHolo():
    def __init__(self, dll_path, verbosity = 2):
        self.dll = ct.cdll.LoadLibrary(os.path.abspath(dll_path))
        self.handleIdx = self.dll.digHoloCreate()
        self.SetVerbosity(verbosity)
        
    def Close(self):
        self.dll.digHoloDestroy(self.handleIdx)
        
    def SetVerbosity(self, verbosity: int):
        """
        Specifies how much information is printed to console/file during operation.
        
        0 : Nothing printed except deep errors from underlying libraries (e.g. MKL). 
        1 : Basic info like errors, start/stop summaries. 
        2 : Debug-like printouts with relatively high level of detail. 
        3 : May God have mercy on your soul.
        Setting non-existent verbosity levels will not raise errors. 
        Specifying less than 0 will default to 0. 
        Specifying above the maximum level will behave like the maximum level.
        
        Parameters
        ----------
        verbosity : int, default 2
            level of verbosity
        """
        ret = self.dll.digHoloConfigSetVerbosity(self.handleIdx, verbosity)
        check_error(ret)
        
    def RedirectConsole(self, file_path = None):
        """
        Calling this function will redirect stdout console to the specified filename.
        
        Parameters
        ----------
        file_path : str or None
            path of the desired output text file.
            If None, redirect the console to stdout
        
        """
        
        self.dll.digHoloConsoleRedirectToFile.argtypes = [ct.c_char_p]
        if file_path:
            charPtr = ct.c_char_p(file_path.encode('utf-8'))
            ret = self.dll.digHoloConsoleRedirectToFile(charPtr)
        else:
            ret = self.dll.digHoloConsoleRestore()
        check_error(ret)
        
    def FrameSimulatorCreateSimple(
        self, 
        frameCount: int, 
        frameWidth: int, 
        frameHeight: int, 
        pixelSize: float, 
        wavelength: float,
        polCount: int = 0, 
        printToConsole: bool = True,
        dataType: str = 'Python'
    ) -> Union[np.array, ct.POINTER(ct.c_float)]:
        """
        A simplified version of the digHoloFrameSimulatorCreate() routine. 
        Whereby the user specifies the minimum necessary parameters to generate frames, and everything else is default.

        This function is also less forgiving than the digHoloFrameSimulatorCreate() routine. e.g. invalid parameters will cause the routine to halt and return a null pointer, rather than choosing default values.
        This function also takes values directly as arguments, rather than pointers. 

        Parameters
        ----------
        frameCount : int
            The number of frames to generate.
        frameWidth : int
            The width of each frame in pixels.
        frameHeight : int
            The height of each frame in pixels.
        pixelSize : float
            The physical dimension of each pixel.
        wavelength: float
             The operating wavelength.
        polCount : int
             The number of polarisation components.
        printToConsole : bool
             Whether or not to print information to the console.
        dataType : {'Python', 'C'}
            Type of data provided, a Python array or a C type array.
            
        Returns
        -------
            frames : ct.POINTER(ct.c_float) or numpy.array
                Simulated frames.
        """
        assert(dataType in ['Python','C'])

        self.dll.digHoloFrameSimulatorCreateSimple.argtypes = [ct.c_int, 
                                                               ct.c_int, 
                                                               ct.c_int,
                                                               ct.c_float,
                                                               ct.c_int,
                                                               ct.c_float,
                                                               ct.c_int]
        self.dll.digHoloFrameSimulatorCreateSimple.restype = ct.POINTER(ct.c_float)
        frameBufferPtr = self.dll.digHoloFrameSimulatorCreateSimple(
            ct.c_int(frameCount),
            ct.c_int(frameWidth),
            ct.c_int(frameHeight),
            ct.c_float(pixelSize),
            ct.c_int(polCount),
            ct.c_float(wavelength),
            ct.c_int(printToConsole)
        )
        if dataType == 'Python':
            return np.ctypeslib.as_array(frameBufferPtr, shape=(frameCount,frameHeight,frameWidth))
        else:
            return frameBufferPtr
        
    def SetFramePixelSize(self, pixelSize: float):
        """
        Sets the current pixel width/height/pitch  of the frame pixels.
        
        Parameters
        ----------
        pixelSize : float
            width=height=pitch of frame pixel
        """
        self.dll.digHoloConfigSetFramePixelSize.argtypes = [ct.c_int, ct.c_float]
        ret = self.dll.digHoloConfigSetFramePixelSize(self.handleIdx, ct.c_float(pixelSize))
        check_error(ret)
        
    def GetFramePixelSize(self) -> float:
        """
        Returns the current pixel width/height/pitch of the frame pixels.
        
        Returns
        -------
        pixelSize : float
            width=height=pitch of frame pixel. 
            Will return 0 for invalid handle index.
        
        """
        self.dll.digHoloConfigGetFramePixelSize.argtypes = [ct.c_int]
        self.dll.digHoloConfigGetFramePixelSize.restype = ct.c_float
        pixelSize = float(self.dll.digHoloConfigGetFramePixelSize(self.handleIdx))
        return pixelSize
        
    def SetFrameDimensions(self, frameWidth: int ,frameHeight: int):
        """
        Sets the full dimensions of the frames in the frame buffer. (e.g. width,height = 640x512)
        
        Parameters
        ----------
        width : int
            Width of the frames. Typically width is the longer dimension on a regular camera. 
            The x-axis. Must be multiple of DIGHOLO_PIXEL_QUANTA.
        height : int
            Height of the frames. Typically the height is the shorter dimension on a regular camera. 
            The y-axis. Must be multiple of DIGHOLO_PIXEL_QUANTA.
        """
        self.dll.digHoloConfigSetFrameDimensions.argtypes = [ct.c_int, ct.c_int];
        ret = self.dll.digHoloConfigSetFrameDimensions(self.handleIdx, frameWidth, frameHeight);
        check_error(ret)
    
    def SetWavelengthCentre(self, lambda0: float):
        """
        Configuring the current operating wavelength. (e.g. for converting angles to k-space)
        Sets the default centre wavelength. 
        The dimensions in Fourier space in terms of angles will depend on the operating wavelength.
        
        Parameters
        ----------
        lambda0 : float
             Operating wavelength.
        """
        self.dll.digHoloConfigSetWavelengthCentre.argtypes = [ct.c_int, ct.c_float]
        ret = self.dll.digHoloConfigSetWavelengthCentre(self.handleIdx, ct.c_float(lambda0))
        check_error(ret)
        
    def GetWavelengthCentre(self) -> float:
        """
        Gets the default centre wavelength. 
        The dimensions in Fourier space in terms of angles will depend on the operating wavelength.
        
        Returns
        -------
        lambda0 : float
             Operating wavelength. Returns zero for invalid handle index.
        """
        self.dll.digHoloConfigGetWavelengthCentre.argtypes = [ct.c_int]
        self.dll.digHoloConfigGetWavelengthCentre.restype = ct.c_float
        lambda0 = float(self.dll.digHoloConfigGetWavelengthCentre(self.handleIdx))
        return lambda0
        
    def SetPolCount(self, polCount: int):
        """
        Sets the number of polarisation components per frame (1 or 2)
        
        Parameters
        ----------
        polCount : float
            Polarisation components per frame (1 or 2).
        """
        self.polCount = polCount
        self.dll.digHoloConfigSetPolCount.argtypes = [ct.c_int]
        ret = self.dll.digHoloConfigSetPolCount(self.handleIdx, polCount)
        check_error(ret)
    
    def SetfftWindowSizeX(self, nx: int):
        ret = self.dll.digHoloConfigSetfftWindowSizeX(self.handleIdx, nx)
        check_error(ret)
        
    def SetfftWindowSizeY(self, ny: int):
        ret = self.dll.digHoloConfigSetfftWindowSizeY(self.handleIdx, ny)
        check_error(ret)
        
    def SetfftWindowSize(self, nx: int, ny: int):
        """
        Sets the width (x-axis) and height (y-axis) of the window within the full-frame which will be FFT'd. 
        Must be multiple of 16 (DIGHOLO_PIXEL_QUANTA).
        For example, the full camera frames may be 640x512, but there is only a 256x256 region
        of interest (window) within that full-frame.
        
        Parameters
        ----------
        width : int
            Width (x-axis) of the FFT window.
        height : int
            Height (y-axis) of the FFT window.
        """
        self.dll.digHoloConfigSetfftWindowSize.argtypes = [ct.c_int, ct.c_int, ct.c_int]
        ret = self.dll.digHoloConfigSetfftWindowSize(self.handleIdx, nx, ny)
        check_error(ret)
        
    def SetIFFTResolutionMode(self, resolutionMode: int):
        ret = self.dll.digHoloConfigSetIFFTResolutionMode(self.handleIdx,resolutionMode)
        check_error(ret)
        
    def SetBasisGroupCount(self, maxMG):
        #Specifies the number of HG mode groups to decompose the beams with
        self.dll.digHoloConfigSetBasisGroupCount(self.handleIdx, maxMG)
        
        
    def SaveConfig(self, file_path: str, polarization: int = 0):
        pixelSize = self.GetFramePixelSize()
        lambda0 = self.GetWavelengthCentre()
        fourierWindowRadius = self.GetFourierWindowRadius()
        beamCentre = self.GetBeamCentre(polarization = polarization)
        defocus = self.GetDefocus(polarization = polarization)
        tilt = self.GetTilt(polarization = polarization)
        np.savez(
            file_path,
            pixelSize = pixelSize,
            lambda0 = lambda0,
            fourierWindowRadius = fourierWindowRadius,
            beamCentre = beamCentre,
            defocus = defocus,
            tilt = tilt,
        )
        
    def LoadConfig(self, file_path: str, polarization: int = 0):
        data = np.load(file_path)
        pixelSize = data['pixelSize']
        lambda0 = data['lambda0']
        fourierWindowRadius = data['fourierWindowRadius']
        beamCentre = data['beamCentre']
        defocus = data['defocus']
        tilt = data['tilt']
        
        self.SetFramePixelSize(pixelSize)
        self.SetWavelengthCentre(lambda0)
        self.SetFourierWindowRadius(fourierWindowRadius)
        self.SetBeamCentre(beamCentre)
        self.SetDefocus(defocus)
        self.SetTilt(tilt)
        
        
    def ConfigOffAxis(self, 
                      frameRes: List[int],
                      fftWindowRes: List[int], 
                      resolutionMode: int, 
                      pixelSize: Union[float, None] = None,
                      lambda0: Union[float, None] = None,
                      maxMG: int = 1,
                      polCount:int = 1
                     ):
        if pixelSize is not None:
            self.SetFramePixelSize(pixelSize)
        if lambda0:
            self.SetWavelengthCentre(lambda0)
        self.SetFrameDimensions(frameRes[0], frameRes[1])
        self.SetfftWindowSize(fftWindowRes[0], fftWindowRes[1])
        self.SetIFFTResolutionMode(resolutionMode)
        self.SetBasisGroupCount(maxMG)
        self.SetPolCount(polCount)
        
    def ConfigSetAutoAlign(
        self,
        enable_align_beam_centre = True,
        enable_align_defocus = True,
        enable_align_tilt = True,
        enable_align_basis_waist = True,
        enable_align_fourier_win_radius = True
        ):
        self.dll.digHoloConfigSetAutoAlignBeamCentre(self.handleIdx, ct.c_int(enable_align_beam_centre))
        self.dll.digHoloConfigSetAutoAlignDefocus(self.handleIdx, ct.c_int(enable_align_defocus))
        self.dll.digHoloConfigSetAutoAlignTilt(self.handleIdx, ct.c_int(enable_align_tilt))
        self.dll.digHoloConfigSetAutoAlignBasisWaist(self.handleIdx, ct.c_int(enable_align_basis_waist))
        self.dll.digHoloConfigSetAutoAlignFourierWindowRadius(self.handleIdx,  ct.c_int(enable_align_fourier_win_radius))
        
    def Convert2Ctypes(self, frames, dtype = ct.c_float):
        frameBufferPtr = np.ctypeslib.as_ctypes(frames.copy().ravel().astype(dtype))
        return frameBufferPtr
    
    def digHoloConfigSetfftWindowSize(self, shape: List[int]) -> None:
        """
        Sets the width (x-axis) and height (y-axis) of the window within the full-frame which will
        be FFT'd. Must be multiple of 16 (DIGHOLO_PIXEL_QUANTA).
        For example, the full camera frames may be 640x512, but there is only a 256x256 region
        of interest (window) within that full-frame.
        """
        width = shape[0]
        height = shape[1]
        
        self.dll.digHoloConfigSetfftWindowSize.argtypes = [ct.c_int, ct.c_int, ct.c_int]
        self.dll.digHoloConfigSetfftWindowSize(self.handleIdx, ct.c_int(width), ct.c_int(height))
        
    def GetfftWindowSize(self):
        """
        Gets the width (x-axis) and height (y-axis) of the window within the full-frame which
        will be FFT'd. Will be multiple of 16 (DIGHOLO_PIXEL_QUANTA).
        For example, the full camera frames may be 640x512, but there is only a 256x256 region
        of interest (window) within that full-frame.
        """
        w = ((ct.c_int))()
        h = ((ct.c_int))()
        
        self.dll.digHoloConfigGetfftWindowSize.argtypes = [
            ct.c_int,
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_int)
        ]
            
        self.dll.digHoloConfigGetfftWindowSize(
            self.handleIdx,
            ct.byref(w),
            ct.byref(h)
        )
        
        return np.int32(w), np.int32(h)
    
    def SetfftWindowSize(self, width, height):
        """
        Sets the width (x-axis) and height (y-axis) of the window within the full-frame which will
        be FFT'd. Must be multiple of 16 (DIGHOLO_PIXEL_QUANTA).
        For example, the full camera frames may be 640x512, but there is only a 256x256 region
        of interest (window) within that full-frame.
        """
        self.dll.digHoloConfigSetfftWindowSize.argtypes = [ct.c_int, ct.c_int, ct.c_int] 
        self.dll.digHoloConfigSetfftWindowSize(self.handleIdx, ct.c_int(width), ct.c_int(height))
    
    def GetFourierWindowRadius(self) -> float:
        """
        Returns the current Fourier window radius in degrees.
        """
        self.dll.digHoloConfigGetFourierWindowRadius.argtypes = [ct.c_int]
        self.dll.digHoloConfigGetFourierWindowRadius.restype = ct.c_float
        return np.float64(self.dll.digHoloConfigGetFourierWindowRadius(self.handleIdx))
    
    
    def SetFourierWindowRadius(self, radius: float):
        """
        Sets the radius (in degrees) of the circular window in Fourier space which selects the
        off-axis term, that will then be IFFT'd to reconstruct the field.
        This is the circular window in Fourier space that is selected. e.g. w_c
        https://doi.org/10.1364/OE.25.033400 Fig. 1 The radius would typically be no more than
        ~1/3 [sqrt(2)/(3+sqrt(2))]=0.320513 of the maximum resolvable angle set by the frame
        pixel size. (e.g. Fig 1c) e.g. if the wavelength is 1565e-9 and the pixelSize is 20e-6,
        w_max would be (1565e-9/(2*20e-6))*(180/pi) = 2.24 degrees, and window radius (w_c)
        should be less than 0.3205*2.24=0.719 degrees The reference beam tilt angle in x and y
        would have to be (x,y)=(3w_c,3w_c)/sqrt(2)=(1.525,1.525) degrees. If the full resolution
        of the camera is not required, smaller windows and shallower reference angles can be
        employed. Smaller windows are also less likely to capture unwanted noise. If Fourier
        wrapping (e.g. Fig 1b) is employed, a larger fractional window could be used [+ root of
        0=8w^2+2-2]->[(-2+sqrt(68))/16]=0.3904. Tilt (x,y) = (w_max-w_c,w_max) Formerly
        called 'digHoloSetApertureSize'
        """ 
        self.dll.digHoloConfigSetFourierWindowRadius.argtypes = [ct.c_int, ct.c_float]
        self.dll.digHoloConfigSetFourierWindowRadius(self.handleIdx, ct.c_float(radius))
    
    def GetIFFTResolutionMode(self):
        """
        Gets whether the IFFT has the same dimensions as the FFT (0), or of the smaller Fourier
        window (1)
        0 : IFFT has the same pixel number and dimensions as the FFT (fftWindowSizeX/Y), and
        hence the reconstructed field is the same dimensions as the original camera frame. This
        mode yields no extra information and is slower, but it is sometimes conveinient. 1 :
        Recommended mode. IFFT has the minimum dimensions, those of the Fourier window.
        Faster with no loss of information.
        """
        self.dll.digHoloConfigGetIFFTResolutionMode.argtypes = [ct.c_int]
        return np.int32(self.dll.digHoloConfigGetIFFTResolutionMode(self.handleIdx))
        
    def SetTilt(self, tilt: List[int], polarization: int = 0):
        """
        Sets the reference beam tilt in degrees. This corresponds with the position of the off-axis term in
        Fourier space that contains our field to be reconstructed.
        """
        self.dll.digHoloConfigSetTilt.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.c_float]
        
        self.dll.digHoloConfigSetTilt(self.handleIdx, ct.c_int(0), ct.c_int(polarization), ct.c_float(tilt[0]))
        self.dll.digHoloConfigSetTilt(self.handleIdx, ct.c_int(1), ct.c_int(polarization), ct.c_float(tilt[1]))
        
    def GetTilt (self, polarization: int = 0) -> List[float]:
        """
        Gets the reference beam tilt in degrees. This corresponds with the position of the off-axis term in
        Fourier space that contains our field to be reconstructed.
        """
        self.dll.digHoloConfigGetTilt.argtypes = [ct.c_int, ct.c_int, ct.c_int]
        self.dll.digHoloConfigGetTilt.restype = ct.c_float
        tilt_x = self.dll.digHoloConfigGetTilt(self.handleIdx, ct.c_int(0), ct.c_int(polarization))
        tilt_y = self.dll.digHoloConfigGetTilt(self.handleIdx, ct.c_int(1), ct.c_int(polarization))
        return [tilt_x, tilt_y]
        
        
    def SetDefocus (self, defocus: float, polarization: int = 0):
        """
        Sets the the reference beam defocus in dioptre.
        """
        self.dll.digHoloConfigSetDefocus.argtypes = [ct.c_int, ct.c_int, ct.c_float]
        self.dll.digHoloConfigSetDefocus (self.handleIdx, ct.c_int(polarization), ct.c_float(defocus))
        
    def GetDefocus(self, polarization: int = 0) -> float:
        """
        Gets the reference beam defocus in dioptre.
        """
        self.dll.digHoloConfigGetDefocus.argtypes = [ct.c_int, ct.c_int]
        self.dll.digHoloConfigGetDefocus.restype = ct.c_float
        return np.float(self.dll.digHoloConfigGetDefocus (self.handleIdx, ct.c_int(polarization)))
    
    def SetBeamCentre(self, center: List[float], polarization: int = 0):
        """
        Sets the current beam centre in the camera plane for a specified axis (x/y) and polarisation
        component.
        """
        self.dll.digHoloConfigSetBeamCentre.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.c_float]
        self.dll.digHoloConfigSetBeamCentre.restype = ct.c_float
        ret = self.dll.digHoloConfigSetBeamCentre(
            self.handleIdx, 
            ct.c_int(0), 
            ct.c_int(polarization), 
            ct.c_float(center[0])
        )
        check_error(ret)
        ret = self.dll.digHoloConfigSetBeamCentre(
            self.handleIdx, 
            ct.c_int(1), 
            ct.c_int(polarization), 
            ct.c_float(center[1])
        )
        check_error(ret)
        
    def GetBeamCentre(self, polarization: int = 0) -> List[float]:
        """
        Returns the current beam centre in the camera plane for a specified axis (x/y) and polarisation
        component.
        """
        self.dll.digHoloConfigGetBeamCentre.argtypes = [ct.c_int, ct.c_int, ct.c_int]
        self.dll.digHoloConfigGetBeamCentre.restype = ct.c_float
        center_x = float(self.dll.digHoloConfigGetBeamCentre(self.handleIdx, ct.c_int(0), ct.c_int(polarization)))
        center_y = float(self.dll.digHoloConfigGetBeamCentre(self.handleIdx, ct.c_int(1), ct.c_int(polarization)))
        return [center_x, center_y]
        
        
    def SetBatch(self, batchCount, frames, dataType = 'Python'):
        """
        Setup a batch of batchCount frames, starting at the frameBufferPtr
        """
        self.batchCount = batchCount
        assert(dataType in ['Python', 'C'])
        if dataType == 'Python':
            self.frameBufferPtr = self.Convert2Ctypes(frames)
            # Needs to be store it otherwise it is freed and creates an error.
            # The C type is created in Python, so Python destroys it
            # when leaving the namespace of the function,  
            # the C array is then destroyed while the dll still needs it.
            # Needs fixing
        else:
            self.frameBufferPtr = frames
        print(self.batchCount)
        self.dll.digHoloSetBatch.argtypes = [ct.c_int, ct.c_int, ct.c_void_p]
        ret = self.dll.digHoloSetBatch(self.handleIdx, batchCount, self.frameBufferPtr)
        check_error(ret)

        
    def ProcessBatch(self):
        """
        Processes a batch of frames using the current settings, 
        and returns a modal decomposition of the resulting fields.
        """
        batchCount = ((ct.c_int))()
        polCount = ((ct.c_int))()
        modeCount = ((ct.c_int))()
            
        self.dll.digHoloProcessBatch.argtypes = [ct.c_int, 
                                                 ct.POINTER(ct.c_int), 
                                                 ct.POINTER(ct.c_int),
                                                 ct.POINTER(ct.c_int)]
        
        ptrOut = self.dll.digHoloProcessBatch(
            self.handleIdx,
            ct.byref(batchCount),
            ct.byref(modeCount),
            ct.byref(polCount)
        )
        batchCount = np.int32(batchCount)
        modeCount = np.int32(modeCount)
        polCount = np.int32(polCount)
        fields = np.ctypeslib.as_array(ptrOut, shape = (batchCount, polCount*modeCount*2))
#         fields = fields[:,0::2]+1j*fields[:,1::2]
        return fields
        
    def AutoAlign(self):
        """
        Run the AutoAlign routine to find parameters like beam centre, tilt, focus and waist.       
        """
        self.dll.digHoloAutoAlign(self.handleIdx)

    def SetRefCalibrationIntensity(self, intensityReference: np.array, wavelengthCount: int = 1):
        """
        Define the intensity of the Reference wave, which will be used to calibrate out a non-uniform
        Reference wave.
        
        Parameters
        ----------
        intensityReference : np.array
            An array of pixels (wavelengthCount x width x height) or (width x height)  specifying
            the intensity of the Reference wave.
        wavelengthCount: int
            The number of frames in 'cal' which is assumed to specify
            Reference waves captured at different wavelengths matching those
            of the FrameBuffer.   
        """
        self.dll.digHoloConfigSetRefCalibrationIntensity.argtypes = [ct.c_int, 
                                                                     ct.c_void_p,
                                                                     ct.c_int,
                                                                     ct.c_int,
                                                                     ct.c_int]
        
        
        width = ct.c_int(intensityReference.shape[-2])
        height = ct.c_int(intensityReference.shape[-1])
        self.referencePtr = self.Convert2Ctypes(intensityReference, dtype = ct.c_ushort)
        ret = self.dll.digHoloConfigSetRefCalibrationIntensity(
            self.handleIdx,
            self.referencePtr,
            ct.c_int(wavelengthCount),
            width,
            height
        )
        check_error(ret)
        
    def SetRefCalibrationField(self, fieldReference: np.array, wavelengthCount : int = 1):
        """
        Define the field of the Reference wave, which will be used to calibrate out a non-uniform and
        aberrated Reference wave.
        See 'digHoloConfigSetRefCalibrationIntensity' for further details. This routine is similar
        to digHoloSetRefCalibrationIntensity, except the user specifies the Reference wave as a
        full field, and hence can compensate not only for non-uniform intensity in the Reference
        wave, but also aberrations.
        
        Parameters
        ----------
        refField : np.array
            An array of pixels (wavelengthCount x width x height) specifying
            the field of the Reference wave
        wavelengthCount : int
            The number of frames in 'cal' which is assumed to specify
            Reference waves captured at different wavelengths matching those
            of the FrameBuffer.
        """
        self.dll.digHoloConfigSetRefCalibrationField.argtypes = [ct.c_int, 
                                                                 ct.c_void_p, 
                                                                 ct.c_int, 
                                                                 ct.c_int,
                                                                 ct.c_int]
        
        width = ct.c_int(fieldReference.shape[-2])
        height = ct.c_int(fieldReference.shape[-1])
        fieldReference = fieldReference.astype(np.complex)
        shape = list(fieldReference.shape)
        shape[-1] *=2
        fieldReferenceFloat = np.zeros(shape, dtype = float)
        fieldReferenceFloat[...,::2] = fieldReference.real
        fieldReferenceFloat[...,1::2] = fieldReference.imag
        self.referenceFieldPtr = self.Convert2Ctypes(fieldReferenceFloat, dtype = ct.c_float)
        
        ret = self.dll.digHoloConfigSetRefCalibrationField(
            self.handleIdx,
            self.referenceFieldPtr,
            wavelengthCount,
            width,
            height
        )
            
    def GetRefCalibrationFields(self):
        """
        Returns Reference wave calibration and corresponding x and y axis.
        
        Returns
        -------
        refFields : np.array
            wavelengthCount x polCount x width x height complex64 array containing the fields
        """
        
        return self.GetRefCalibrationFieldsAndAxis()[0]
    
    def GetRefCalibrationFieldsAndAxis(self):
        """
        Returns Reference wave calibration and corresponding x and y axis.
        
        Returns
        -------
        refFields : np.array
            wavelengthCount x polCount x width x height complex64 array containing the fields
        X : np.array
            x-axis of the field.
        Y : np.array
            y-axis of the field.
        
        """
        wavelengthCount = ((ct.c_int))()
        polCount = ((ct.c_int))()
        
        #The x/y axis of the field. Corresponding with the dimension in the camera
        #plane.
        xPtr = (ct.POINTER(ct.c_float))()
        yPtr = (ct.POINTER(ct.c_float))()
        #The width/height of the x and y axes respectively.
        width = ((ct.c_int))()
        height = ((ct.c_int))()
        
        self.dll.digHoloConfigGetRefCalibrationFields.argtypes = [ct.c_int, 
                                                                  ct.POINTER(ct.c_int), 
                                                                  ct.POINTER(ct.c_int), 
                                                                  ct.POINTER(ct.POINTER(ct.c_float)),
                                                                  ct.POINTER(ct.POINTER(ct.c_float)),
                                                                  ct.POINTER(ct.c_int),
                                                                  ct.POINTER(ct.c_int)]
        
        self.dll.digHoloConfigGetRefCalibrationFields.restype = ct.POINTER(ct.c_float)
        
        ptrOut = self.dll.digHoloConfigGetRefCalibrationFields(
            self.handleIdx,
            ct.byref(wavelengthCount),
            ct.byref(polCount),
            ct.byref(xPtr),
            ct.byref(yPtr),
            ct.byref(width),
            ct.byref(height)
        )
        
                
        wavelengthCount = np.int32(wavelengthCount)
        #The number of polarisation components per frame (batch element)
        polCount = np.int32(polCount)
        #The width/height of the reconstructed field per polarisation
        width = np.int32(width)
        height = np.int32(height)

        X = np.ctypeslib.as_array(xPtr,shape=[width])
        Y = np.ctypeslib.as_array(yPtr,shape=[height])
        
        refFields = np.ctypeslib.as_array(ptrOut, shape = (wavelengthCount, polCount, width, height*2)).squeeze()
        refFields = refFields[...,0::2]+1j*refFields[...,1::2]
        return refFields, X, Y
                                                                  
    def GetFields(self):
        """
        Return reconstructed fields.
        
        Returns
        -------
        fields : np.array
            The function returns the batchCount (number of fields) and the polCount
            (number of polarisation components per field).
        """
        return self.GetFieldsAndAxis()[0]
        
    def GetFieldsAndAxis(self):
        """
        Return reconstructed fields and corresponding x and y axis.
        
        
        Returns
        -------
        fields : np.array
            The function returns the batchCount (number of fields) and the polCount
            (number of polarisation components per field).
        X : np.array
            x-axis of the field.
        Y : np.array
            y-axis of the field.
        """
        assert(self.batchCount)
        batchCount = ((ct.c_int))()
        polCount = ((ct.c_int))()

        #The x/y axis of the field. Corresponding with the dimension in the camera
        #plane.
        xPtr = (ct.POINTER(ct.c_float))()
        yPtr = (ct.POINTER(ct.c_float))()
        #The width/height of the x and y axes respectively.
        width = ((ct.c_int))()
        height = ((ct.c_int))()

        self.dll.digHoloGetFields.argtypes = [ct.c_int,
                                             ct.POINTER(ct.c_int),
                                             ct.POINTER(ct.c_int), 
                                             ct.POINTER(ct.POINTER(ct.c_float)),
                                             ct.POINTER(ct.POINTER(ct.c_float)),
                                             ct.POINTER(ct.c_int),
                                             ct.POINTER(ct.c_int)]
        #Output return types
        self.dll.digHoloGetFields.restype = ct.POINTER(ct.c_float)

        ptrOut = self.dll.digHoloGetFields(
            self.handleIdx,
            ct.byref(batchCount),
            ct.byref(polCount),
            ct.byref(xPtr),
            ct.byref(yPtr),
            ct.byref(width),
            ct.byref(height)
        )
        #The number of camera frames in the batch returned by digHoloGetFields
        batchCount = np.int32(batchCount)
        #The number of polarisation components per frame (batch element)
        polCount = np.int32(polCount)
        #The width/height of the reconstructed field per polarisation
        width = np.int32(width)
        height = np.int32(height)

        fields = np.ctypeslib.as_array(ptrOut, shape = (batchCount, polCount*width, height*2))
        fields = fields[:,:,0::2]+1j*fields[:,:,1::2]
        X = np.ctypeslib.as_array(xPtr,shape=[width])
        Y = np.ctypeslib.as_array(yPtr,shape=[height])
        return fields, X, Y
