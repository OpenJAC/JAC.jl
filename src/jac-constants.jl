#=
2014 CODATA recommended values, obtained from http://physics.nist.gov/cuu/Constants/
=#

# Dimensionless constants
const FINE_STRUCTURE_CONSTANT         = 7.297_352_566_4e-3 
const INVERSE_FINE_STRUCTURE_CONSTANT = 137.035_999_139

# Constants in SI units
const BOHR_RADIUS_SI                  = 0.529_177_210_67e-10
const BOLTZMANN_CONSTANT_SI           = 1.380_648_52e-23
const ELECTRON_MASS_SI                = 9.109_383_56e-31
const ELEMENTARY_CHARGE_SI            = 1.602_176_620_8e-19
const HARTREE_ENERGY_SI               = 4.359_744_650e-18
const PLANCK_CONSTANT_SI              = 6.626_070_040e-34
const PLANCK_CONSTANT_OVER_2_PI_SI    = 1.054_571_800e-34
const SPEED_OF_LIGHT_IN_VACUUM_SI     = 299_792_458

# Constants with eV instead of J
const BOLTZMANN_CONSTANT_EV           = 8.617_330_3e-5
const HARTREE_ENERGY_EV               = 27.211_386_02
const PLANCK_CONSTANT_EV              = 4.135_667_662e-15
const PLANCK_CONSTANT_OVER_2_PI_EV    = 6.582_119_514e-16

# Constants with other atomic units
const RYDBERG_IN_KAYSERS              = 109.73731534e3
const ELECTRON_MASS_IN_G              = 9.1093897e-28
const HBAR_IN_ERGS                    = 1.05457266e-27
const ELECTRON_CHARGE_IN_ESU          = 4.80320680e-10
   
# Constants in unified atomic mass units (u)
const ELECTRON_MASS_U                 = 5.485_799_090_70e-4
const NEUTRON_MASS_U                  = 1.008_664_915_88
const PROTON_MASS_U                   = 1.007_276_466_879

# Relationships between energy equivalents
const ELECTRON_VOLT_ATOMIC_MASS_UNIT_RELATIONSHIP = 1.073_544_110_5e-9

# Predefined conversion factors
const CONVERT_ENERGY_AU_TO_EV           = HARTREE_ENERGY_EV
const CONVERT_ENERGY_AU_TO_KAYSERS      = 2.0 * RYDBERG_IN_KAYSERS
const CONVERT_ENERGY_AU_TO_PER_SEC      = 4.134138e16
const CONVERT_TIME_AU_TO_SEC            = 2.418_884_254e-17
const CONVERT_CROSS_SECTION_AU_TO_BARN  = BOHR_RADIUS_SI^2 * 10.0e28
const CONVERT_RATE_AU_TO_PER_SEC        = (ELECTRON_MASS_IN_G/HBAR_IN_ERGS) * ((ELECTRON_CHARGE_IN_ESU^2/HBAR_IN_ERGS)^2)
const CONVERT_LENGTH_AU_TO_FEMTOMETER   = BOHR_RADIUS_SI * 1.0e15

# Predefined coefficients for numerical integration
const FINITE_DIFFERENCE_NPOINTS         = 6
const fivePointCoefficients             = Array{Float64}([2 * 7, 32, 12, 32, 7]) * 2 / 45  # * grid.h     
const newtonCotesCoefficients           = fivePointCoefficients

n = JAC.FINITE_DIFFERENCE_NPOINTS     
##x const weights = Array{Float64}(2*n + 1, 2*n + 1)     
const weights = Array{Float64}(undef, 2*n + 1, 2*n + 1)
for i = 1:2*n + 1       
    weights[i,:] = JAC.Math.finiteDifferenceWeights(-n + i - 1, 2*n + 1, order=1)[2,:]     
end


#
#
# Global settings that can be (re-) defined by the user.
JAC_FRAMEWORK               = "relativistic"
JAC_CONT_SOLUTION           = ContBessel                ###  ContBessel    ContSine     AsymptoticCoulomb    NonrelativisticCoulomb    BsplineGalerkin
JAC_CONT_NORMALIZATION      = PureSine                  ###  PureSine    CoulombSine    OngRussek
JAC_QED_HYDROGENIC_LAMBDAC  = [1.0,  1.0,  1.0,  1.0,  1.0]
JAC_QED_NUCLEAR_CHARGE      = 0.1
JAC_QED_MODEL               = QedPetersburg             ###  QedPetersburg  QedSydney    
##x JAC_CONT_PHASE              = WrtSin                ###  WrtSin
JAC_WARNINGS                = JAC_WARNINGS = String[]

JAC_ENERGY_UNIT             = "eV"
JAC_CROSS_SECTION_UNIT      = "barn"
JAC_RATE_UNIT               = "1/s"
JAC_TIME_UNIT               = "sec"

JAC_SUMMARY_IOSTREAM        = stdout
JAC_TEST_IOSTREAM           = stdout
JAC_PRINT_SUMMARY           = false
JAC_PRINT_TEST              = false

# global JAC_STANDARD_SUBSHELL_LIST  = nothing
# JAC_STANDARD_GRID           = Radial.Grid("grid: exponential")

