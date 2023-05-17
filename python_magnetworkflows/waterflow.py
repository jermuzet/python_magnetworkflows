from dataclasses import dataclass

import math
import json

import warnings
from pint import UnitRegistry, Unit, Quantity

# Ignore warning for pint
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])

import freesteam

# Pint configuration
ureg = UnitRegistry()
ureg.default_system = "SI"
ureg.autoconvert_offset_to_baseunit = True


def rho(Tw: float, P: float) -> float:
    """
    compute water volumic mass in ???

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'rho={Steam.rho}')
    return Steam.rho


def Cp(Tw: float, P: float) -> float:
    """
    compute water specific heat in ???

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'Cp={Steam.cp}')
    return Steam.cp


@dataclass
class waterflow:
    Vpump0: float = 1000  # rpm
    Vpmax: float = 2840  # rpm
    F0_l_per_second: float = 0  # l/s
    Fmax_l_per_second: float = 140  # l/s
    Pmax: float = 22  # bar
    Pmin: float = 4  # bar
    Imax: float = 28000  # A

    # Flow params
    @classmethod
    def flow_params(cls, filename: str):
        with open(filename, "r") as f:
            flow_params = json.loads(f.read())

        Vpump0 = flow_params["Vp0"]["value"]  # rpm
        Vpmax = flow_params["Vpmax"]["value"]  # rpm
        F0_l_per_second = flow_params["F0"]["value"]  # l/s
        Fmax_l_per_second = flow_params["Fmax"]["value"]  # l/s
        Pmax = flow_params["Pmax"]["value"]  # bar
        Pmin = flow_params["Pmin"]["value"]  # bar
        Imax = flow_params["Imax"]["value"]  # Amperes

        return cls(Vpump0, Vpmax, F0_l_per_second, Fmax_l_per_second, Pmax, Pmin, Imax)

    def vpump(self, objectif: float) -> float:
        Vpump = self.Vpmax
        if objectif <= self.Imax:
            Vpump = (
                self.Vpump0 + (self.Vpmax - self.Vpump0) * (objectif / self.Imax) ** 2
            )

        return Vpump

    def flow(self, objectif: float) -> float:
        """
        compute flow in m^3/s
        """

        units = [
            ureg.liter / ureg.second,
            ureg.meter * ureg.meter * ureg.meter / ureg.second,
        ]
        F0 = Quantity(self.F0_l_per_second, units[0]).to(units[1]).magnitude
        Fmax = Quantity(self.Fmax_l_per_second, units[0]).to(units[1]).magnitude
        return F0 + Fmax * self.vpump(objectif) / self.Vpmax

    def pressure(self, objectif: float) -> float:
        """
        compute pressure in bar ???
        """
        return (
            self.Pmin
            + (self.Pmax - self.Pmin) * (self.vpump(objectif) / self.Vpmax) ** 2
        )

    def umean(self, objectif: float, section: float) -> float:
        """
        compute umean in m/s ???
        """
        # print("flow:", flow(objectif), section)
        return self.flow(objectif) / section

    def montgomery(self, Tw: float, Umean: float, Dh: float) -> float:
        """
        compute heat exchange coefficient in ??

        Tw: K
        Umean: m/s
        Dh: meter
        """
        return (
            1426
            * (1 + 1.5e-2 * (Tw - 273))
            * math.exp(math.log(Umean) * 0.8)
            / math.exp(math.log(Dh) * 0.2)
        )
