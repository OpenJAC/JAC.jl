
"""
`module  JAC.Basics`  
	... a submodel of JAC that contains many basic types/struct that are specific to the JAC module; this module also defines
	    a number of basic functions/methods that are later extended. We here provide proper docstrings for all abstract and
	    conrete types (struct) in order to allow the user to easy understand and access the individual fields of each type/struct 
	    definition.
"""
module Basics


using Printf, QuadGK

include("module-Basics-inc-first.jl")
include("module-Basics-inc-abstract.jl")
include("module-Basics-inc-second.jl")

export  add, analyze, compute, diagonalize, estimate, generate, interpolate, integrate, modify, perform, provide, tabulate       
        
# Functions/methods that are later added to the module Basics
function add                                                    end
function addZerosToCsfR                                         end
function analyze                                                end
function analyzeConvergence                                     end
function analyzeGrid                                            end
function compute                                                end
function computeDensity                                         end
function computeDiracEnergy                                     end
function computeMeanSubshellOccupation                          end
function computeMultipletForGreenApproach                       end
function computePotential                                       end
function computeScfCoefficients                                 end
function determineEnergySharings                                end
function determineHoleShells                                    end
function determineMeanEnergy                                    end
function determineNearestPoints                                 end
function determineNonorthogonalShellOverlap                     end
function determineParity                                        end
function determinePolarizationLambda                            end
function determinePolarizationVector                            end
function determineSelectedLines                                 end
function determineSelectedPathways                              end
function diagonalize                                            end
function diracDelta                                             end
function display                                                end
function displayLevels                                          end
function displayMeanEnergies                                    end
function displayMeshes                                          end
function displayOpenShells                                      end
function displayOrbitalOverlap                                  end
function displayOrbitalProperties                               end
function excludeConfigurations                                  end
function excludeDoubles                                         end
function expandOrbital                                          end
function extractLeadingConfiguration                            end
function extractLeadingConfigurationR                           end
function extractMeanEnergy                                      end
function extractNoOpenShells                                    end
function extractNonrelativisticShellList                        end
function extractNonrelativisticConfigurations                   end
function extractNonrelativisticConfigurationFromCsfR            end
function extractOpenShells                                      end
function extractOpenSubshells                                   end
function extractOpenShellQNfromCsfNR                            end
function extractOpenShellQNfromCsfR                             end
function extractRelativisticConfigurations                      end
function extractRelativisticConfigurationFromCsfR               end
function extractRelativisticSubshellList                        end
function extractRydbergSubshellList                             end
function extractShellList                                       end
function extractShellOccupationDifference                       end
function extractShellOccupationFromCsfR                         end
function extractSubshellList                                    end
function extractValenceShell                                    end
function FermiDirac                                             end
function generate                                               end
function generateBasis                                          end
function generateConfigurationRs                                end
function generateConfigurations                                 end
function generateConfigurationsForExcitationScheme              end
function generateConfigurationsWithAdditionalElectron           end
function generateConfigurationsWithAdditionalElectrons          end
function generateConfigurationsWithElectronCapture              end
function generateConfigurationsWithElectronLoss                 end
function generateCsfRs                                          end
function generateFieldCoordinates                               end
function generateLevelWithExtraElectron                         end
function generateLevelWithExtraSubshell                         end
function generateLevelWithExtraSubshells                        end
function generateLevelWithExtraTwoElectrons                     end
function generateLevelWithSymmetryReducedBasis                  end
function generateMeshCoordinates                                end
function generateOrbitalsForBasis                               end
function generateOrbitalsForPotential                           end
function generateOrbitalSuperposition                           end
function generateShellList                                      end
function generateSubshellList                                   end
function generateSpectrumLorentzian                             end
function generateSpectrumGaussian                               end
function hasSubshell                                            end
function integrate                                              end
function integrateOnGridNewtonCotes                             end
function integrateOnGridSimpsonRule                             end
function integrateOnGridTrapezRule                              end
function interpolateOnGridGrasp92                               end
function interpolateOnGridTrapezRule                            end
function isSimilar                                              end
function isSymmetric                                            end
function isStandardSubshellList                                 end
function isViolated                                             end
function isZero                                                 end
function lastPoint                                              end
function merge                                                  end
function modifyLevelEnergies                                    end
function modifyLevelMixing                                      end
function perform                                                end
function performSCF                                             end
function performCI                                              end
function plot                                                   end
function read                                                   end
function readCslFileGrasp92                                     end
function readOrbitalFileGrasp92                                 end
function readMixFileRelci                                       end
function readMixingFileGrasp92                                  end
function recast                                                 end
function selectLevel                                            end
function selectLevelPair                                        end
function selectLevelTriple                                      end
function selectSymmetry                                         end
function shiftTotalEnergies                                     end
function sortByEnergy                                           end
function subshellStateString                                    end
function tabulate                                               end
function tabulateKappaSymmetryEnergiesDirac                     end
function tools                                                  end
function yesno                                                  end

end  ## module
