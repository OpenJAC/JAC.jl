
"""
`module  JAC.HighHarmonic
    ... a submodel of JAC that contains all methods to set-up and process high-harmonic computations. 
        It defines in particular a number of structs to deal with the observables, approaches, laser pulses and
        time-dependent fields as they frequently occur in the computation of high-harmonic spectra.
"""
module HighHarmonic

    using Printf, ..Basics, ..Defaults, ..Radial, ..ManyElectron, ..Nuclear, ..TableStrings

    
    """
    `abstract type HighHarmonic.AbstractTargetModel` 
        ... defines an abstract and a number of other types to represent an atomic target (cloud) for HHG:

      + struct LocalizedTargetModel     ... to represent a point-like localized target       
      + struct GaussianTargetModel      ... to represent a localized, Gaussian-type target (not yet)      
      + struct PlanarTargetModel        ... to represent a planar target of finite extent (not yet)      
      + struct UserTargetModel          ... to represent a user-defined 3-dimensional target (not yet)
    """
    abstract type  AbstractTargetModel  end

    
    """
    `struct  HighHarmonic.LocalizedTargetModel  <:  AbstractTargetModel`  
        ... defines a struct for point-like localized target cloud.

        + center     ::CartesianPoint   ... Center of the (point-like) target.
    """
    struct  LocalizedTargetModel  <:  AbstractTargetModel
        center       ::CartesianPoint
    end

    
    """
    `struct  HighHarmonic.GaussianTargetModel  <:  AbstractTargetModel`  
        ... defines a struct for localized, Gaussian target cloud.

        + center     ::CartesianPoint   ... Center of the (point-like) target.
        + NoPoints   ::Int64            ... Number of points to be used for the target.
        + sigma      ::Float64          ... Gaussian half-widths of the cloud.
    """
    struct  GaussianTargetModel  <:  AbstractTargetModel
        center       ::CartesianPoint
        NoPoints     ::Int64  
        sigma        ::Float64 
    end

    
    """
    `struct  HighHarmonic.TargetCloud`  
        ... defines a struct to represent a -- point-like or extended -- target cloud for the computation of high-harmonic spectra

        + model      ::AbstractTargetModel           ... Model of the target cloud, from which the target is generated.
        + points     ::Array{WeightedCartesian,1}    ... Cartesian points and weights that define the atomic target cloud.
    """
    struct  TargetCloud
        model        ::AbstractTargetModel
        points       ::Array{WeightedCartesian,1}
    end

    
    
    """
    `abstract type HighHarmonic.AbstractDeterctorModel` 
        ... defines an abstract and a number of other types to represent an dectector screen for high-harmonic spectra:

      + struct LocalizedDetectorModel      ... to represent a point-like localized detector       
      + struct SphericalDetectorModel      ... to represent a partial-sphere detector screen (not yet)      
      + struct PlanarDetectorModel         ... to represent a planar detector screen of finite extent (not yet)      
    """
    abstract type  AbstractDetectorModel  end

    
    """
    `struct  HighHarmonic.LocalizedDetectorModel  <:  AbstractDetectorModel`  
        ... defines a struct for point-like localized detector.

        + center     ::CartesianPoint   ... Center of the (point-like) detector.
    """
    struct  LocalizedDetectorModel  <:  AbstractDetectorModel
        center       ::CartesianPoint
    end

    
    """
    `struct  HighHarmonic.SphericalDetectorModel  <:  AbstractDetectorModel`  
        ... defines a struct for a partly-spherical detector screen.

        + center     ::CartesianPoint   ... Center of the detector.
        + NoPoints   ::Int64            ... Number of points to be used for the detector screen.
        + theta      ::Float64          ... polar (arc) angle of the screen.
    """
    struct  SphericalDetectorModel  <:  AbstractDetectorModel
        center       ::CartesianPoint
        NoPoints     ::Int64  
        theta        ::Float64 
    end

    
    """
    `struct  HighHarmonic.DetectorScreen`  
        ... defines a struct to represent a -- point-like or extended -- detector screen for the computation of high-harmonic spectra.

        + model      ::AbstractDetectorModel           ... Model of the detector, from which the detection screen is generated.
        + points     ::Array{WeightedCartesian,1}      ... Cartesian points and weights that define the detection screen.
    """
    struct  DetectorScreen
        model        ::AbstractDetectorModel
        points       ::Array{WeightedCartesian,1}
    end
    
    

    """
    `abstract type HighHarmonic.AbstractHhgObservable` 
        ... defines an abstract and a number of singleton types for the observables that can be computed for HHG processes.

      + struct UniColorSpectrum     
        ... to compute a standard (single-color) high-harmonic spectrum (for a point-like target and point-like detector).       
            
      + struct BiColorSpectrum
        ... to compute a bi-color high-harmonic spectrum (for a point-like target and point-like detector).
            
      + struct AttoSecondPulse     (not yet)
      + struct IntensityProfile    (not yet)
      + struct PhaseProfile        (not yet)
            
      + struct UniColorSpectrum
    """
    abstract type  AbstractHhgObservable                         end
    struct         UniColorSpectrum  <:  AbstractHhgObservable   end
    struct         BiColorSpectrum   <:  AbstractHhgObservable   end
    
    

    """
    `abstract type HighHarmonic.AbstractHhgApproach` 
        ... defines an abstract and a number of singleton types for the computational approach/model that is applied in the
            computation of HHG spectra.

      + struct HhgHydrogenicSaddlePoint         
        ... to apply hydrogenic orbitals for the computation of the time-dependent dipole amplitude and the saddle-point
            approximation for the dipole moments.
            
      + struct HhgHydrogenicFullVolkov
        ... to apply hydrogenic orbitals for the computation of the time-dependent dipole amplitude but the full
             Volkov phase for the dipole moments.
    """
    abstract type  AbstractHhgApproach                                 end
    struct         HhgHydrogenicSaddlePoint  <:  AbstractHhgApproach   end
    struct         HhgHydrogenicFullVolkov   <:  AbstractHhgApproach   end

    
    """
    `abstract type HighHarmonic.AbstractHhgPulse` 
        ... defines an abstract and a number of other types to represent an incident laser pulse/field for HHG:

      + struct HhgPlaneWaveLaser        ... to represent a long, plane-wave laser field 
      + struct HhgSineSquaredPulse      ... to represent a sin^2(...) pulse  (not yet)      
    """
    abstract type  AbstractHhgPulse  end

    
    """
    `struct  HighHarmonic.HhgPlaneWaveLaser  <:  AbstractHhgPulse`  
        ... defines a struct for a long, plane-wave laser field of given frequency and intensity.

        + omega          ::Float64   ... Frequency of the incident laser field.
        + intensity      ::Float64   ... Intensity [W/cm^2] of the incident laser field.
        + polarization   ::...       ... Polarization of the incident laser field.
        + beamProfile    ::...
        + envelope       ::...
        + length         ::...
        + centralLength  ::...
    """
    struct  HhgPlaneWaveLaser  <:  AbstractHhgPulse
        omega            ::Float64
        intensity        ::Float64
    end


    """
    `HighHarmonic.HhgPlaneWaveLaser()`  ... constructor for an (empty instance of .
    """
    function HhgPlaneWaveLaser()
        HhgPlaneWaveLaser(0., 0.)
    end


    # `Base.show(io::IO, pulse::HhgPlaneWaveLaser)`  ... prepares a proper printout of the variable pulse::HhgPlaneWaveLaser.
    function Base.show(io::IO, pulse::HhgPlaneWaveLaser) 
        println(io, "omega:                     $(pulse.omega)  ")
        println(io, "intensity:                 $(pulse.intensity)  ")
    end
    


    """
    `struct  HighHarmonic.Computation`  
        ... defines a type for the computation of high-harmonic spectra, atto-second pulses, intensity and phase profiles,
            and several others.

        + observable         ::AbstractHhgObservable          ... Specifies one of the observables to be calculated.
        + approach           ::AbstractHhgApproach            ... Specifies the approach/frame to be applied in the computations.
        + pulse              ::AbstractHhgPulse               ... Specifies the incident laser field.   
        + nuclearModel       ::Nuclear.Model                  ... Model, charge and parameters of the nucleus.
        + initialOrbital     ::Radial.Orbital                 ... (Single-electron) orbital to specify the initial state of atom.
        + target             ::TargetCloud                    ... Specifies the target for generating high harmonics.
        + detector           ::DetectorScreen                 ... Specifies the detector for observing high harmonics.
    """
    struct  Computation
        observable           ::AbstractHhgObservable
        approach             ::AbstractHhgApproach
        pulse                ::AbstractHhgPulse
        nuclearModel         ::Nuclear.Model
        initialOrbital       ::Radial.Orbital
        target               ::TargetCloud 
        detector             ::DetectorScreen 
    end 


    """
    `HighHarmonic.Computation()`  ... constructor for an 'empty' instance of a HighHarmonic.Computation.
    """
    function Computation()
        Computation("",  Nuclear.Model(0.) )
    end


    # `Base.show(io::IO, computation::HighHarmonic.Computation)`  ... prepares a proper printout of the variable computation::HighHarmonic.Computation.
    function Base.show(io::IO, computation::HighHarmonic.Computation) 
        println(io, "name:                     $(computation.name)  ")
    end


    """
    `HighHarmonic.computeFrequencyDipoleMoment()` 
        ... to compute the frequency-dependent dipole moment D^~ (omega)
    """
    function computeFrequencyDipoleMoment()
        return( 0. )
    end


    """
    `HighHarmonic.computeTimeDipoleAmplitude()` 
        ... to compute the time-dependent dipole amplitude d(t)
    """
    function computeTimeDipoleAmplitude()
        return( 0. )
    end


    """
    `HighHarmonic.computeTimeDipoleMoment()` 
        ... to compute the time-dependent dipole moment D(t)
    """
    function computeTimeDipoleMoment()
        return( 0. )
    end


    """
    `HighHarmonic.computeFieldAndPotential()` 
        ... to compute the electric field and the vector potential for the given pulse at time t
    """
    function computeFieldAndPotential()
        return( 0. )
    end


    """
    `HighHarmonic.perform(comp::HighHarmonic.Computation)` 
        ... to perform the computation
    """
    function perform(comp::HighHarmonic.Computation)
        return( 0. )
    end

end # module


