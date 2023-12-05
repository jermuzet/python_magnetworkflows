from dataclasses import dataclass

import json

import warnings
from pint import UnitRegistry, Unit, Quantity

# Ignore warning for pint
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


# Pint configuration
ureg = UnitRegistry()
ureg.default_system = "SI"
ureg.autoconvert_offset_to_baseunit = True


@dataclass
class waterflow:
    Vpump0: float = 1000  # rpm
    Vpmax: float = 2840  # rpm
    F0_l_per_second: float = 0  # l/s
    Fmax_l_per_second: float = 140  # l/s
    Pmax: float = 22  # bar
    Pmin: float = 4  # bar
    BP: float = 4  # bar
    Imax: float = 28000  # A

    # Flow params
    @classmethod
    def flow_params(cls, pwd: str, filename: str):
        with open(f"{pwd}/{filename}", "r") as f:
            flow_params = json.loads(f.read())

        Vpump0 = flow_params["Vp0"]["value"]  # rpm
        Vpmax = flow_params["Vpmax"]["value"]  # rpm
        F0_l_per_second = flow_params["F0"]["value"]  # l/s
        Fmax_l_per_second = flow_params["Fmax"]["value"]  # l/s
        Pmax = flow_params["Pmax"]["value"]  # bar
        Pmin = flow_params["Pmin"]["value"]  # bar
        BP = 4
        if "BP" in flow_params:
            BP = flow_params["BP"]["value"]  # bar
        Imax = flow_params["Imax"]["value"]  # Amperes

        return cls(
            Vpump0, Vpmax, F0_l_per_second, Fmax_l_per_second, Pmax, Pmin, BP, Imax
        )

    def vpump(self, objectif: float) -> float:
        Vpump = self.Vpmax + self.Vpump0
        if objectif <= self.Imax:
            Vpump = self.Vpmax * (objectif / self.Imax) ** 2 + self.Vpump0

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
        F = F0 + Fmax
        if objectif <= self.Imax:
            F = F0 + Fmax * self.vpump(objectif) / (self.Vpmax + self.Vpump0)
        return F

    def pressure(self, objectif: float) -> float:
        """
        compute pressure in bar
        """

        # print(f'pressure: P(I={objectif})={self.pressure(objectif)}, P0={self.Pmin}, Pmax={self.Pmax}')

        if objectif <= self.Imax:
            return (
                self.Pmin
                + self.Pmax * (self.vpump(objectif) / (self.Vpmax + self.Vpump0)) ** 2
            )
        return self.Pmin + self.Pmax

    def dpressure(self, objectif: float) -> float:
        """
        compute dpressure in bar ???
        """
        # print(f'dpressure: P(I={objectif})={self.pressure(objectif)}, BP={self.BP}')
        return self.pressure(objectif) - self.BP

    def umean(self, objectif: float, section: float) -> float:
        """
        compute umean in m/s ???
        """
        # print("flow:", flow(objectif), section)
        return self.flow(objectif) / section
