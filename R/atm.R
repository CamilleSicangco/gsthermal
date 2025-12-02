#' Calculate net radiation
#'
#' @param Tair Air temperature, deg C
#' @param PPFD Photosynthetic photon flux density, mu mol m-2 s-1
#' @param VPD Atmospheric vapor pressure deficit, kPa
#' @param albedo Leaf albedo, unitless
#' @param epsilon Leaf emissivity, unitless
#'
#' @return Net radiation, W m-2
#' @export
calc_Rnet = function(Tair = 25,
                     PPFD = 1000,
                     VPD = 1.5,
                     albedo = 0.15,
                     epsilon = 0.97
)
{
  # Constants
  SW_2_PAR = 4.57 * 0.5  # SW (W m-2) to PAR (mumol m-2 s-1)
  PAR_2_SW = 1 / SW_2_PAR
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m-2 K-4

  # unit conversions
  TairK = Tair + 273 # air temperature in Kelvin
  emissivity = calc_emissivity(Tair, VPD)
  # incoming short and long wave radiation
  Rsw = (1. - albedo) * PPFD * PAR_2_SW  # W m-2
  Rlw = emissivity * sigma * TairK ** 4.  # W m-2

  Rnet = Rsw + Rlw - epsilon * sigma * TairK ** 4 # W m-2
  return(Rnet)
}


#' Calculate atmospheric emissivity
#'
#' @param Tair Air temperature, deg C
#' @param VPD Atmospheric vapor pressure deficit, kPa
#'
#' @return Apparent atmospheric emissivity, unitless
#' @noRd
calc_emissivity = function(Tair = 25,
                           VPD = 1.5
)
{
  # Constants
  sigma = 5.67e-8 # Stefan-Boltzman constant, W m-2 K-4

  TairK = Tair + 273 # air temperature in Kelvin
  ea = (plantecophys::esat(Tair) - 4)*1e3
  emissivity = (0.031 * ea + 2.84 * TairK - 522.5) / (sigma * TairK ** 4.)
  return(emissivity)
}
