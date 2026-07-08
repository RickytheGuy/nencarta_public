try:
    from enum import StrEnum
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):
        pass

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

class NencartaEnum(StrEnum):
    """
    Base class for all nencarta enumerations. Inherit from this instead of StrEnum to ensure that the enum values are always strings, which is important for serialization and user-friendliness.
    """
    @classmethod
    def from_string(cls, value: str) -> Self:
        try:
            return cls(value.lower())
        except ValueError:
            raise ValueError(
                f"Invalid value '{value}'. "
                f"Valid values: {[m.value for m in cls]}"
            )
        
    @classmethod
    def from_string_exact_case(cls, value: str) -> Self:
        try:
            return cls(value)
        except ValueError:
            raise ValueError(
                f"Invalid value '{value}'. "
                f"Valid values: {[m.value for m in cls]}"
            )

class FloodMapMode(NencartaEnum):
    FORECAST = 'forecast'
    USER = 'user'
    RETURN_PERIOD = 'return_period'

class StreamflowSource(NencartaEnum):
    GEOGLOWS = 'geoglows'
    NWM_SHORT_RANGE = 'nwm_short_range'
    NWM_MEDIUM_RANGE = 'nwm_medium_range'
    NWM_LONG_RANGE = 'nwm_long_range'

    def is_nwm(self) -> bool:
        return self in {
            self.NWM_SHORT_RANGE,
            self.NWM_MEDIUM_RANGE,
            self.NWM_LONG_RANGE,
        }

    def __str__(self):
        if self == self.GEOGLOWS:
            return self.value.upper()
        return self.value[:3].upper() + self.value[3:]
    
    @property
    def upstream_id(self) -> str:
        if self == self.GEOGLOWS:
            return "LINKNO"
        return "COMID"

    @property
    def downstream_id(self) -> str:
        if self == self.GEOGLOWS:
            return "DSLINKNO"
        return "TOCOMID"

class Mapper(NencartaEnum):
    FLOODSPREADER = 'FloodSpreader'
    CURVE2FLOOD_KERNEL_WEIGHTED = 'Curve2Flood-Kernel Weighted'
    CURVE2FLOOD_FLDPLNPY = 'Curve2Flood-FLDPLNpy'
    CURVE2FLOOD_MULTI_POINT_INTERPOLATION = 'Curve2Flood-Multi-Point Interpolation'

    @classmethod
    def from_string(cls, value):
        return cls.from_string_exact_case(value)
    
    @classmethod
    def from_string_exact_case(cls, value):
        if value is None:
            return Mapper.CURVE2FLOOD_KERNEL_WEIGHTED

        value = str(value).strip()
        if value == "":
            return Mapper.CURVE2FLOOD_KERNEL_WEIGHTED
        
        return super().from_string_exact_case(value)
    
    def is_curve2flood_mapper(self) -> bool:
        return self in {
            self.CURVE2FLOOD_KERNEL_WEIGHTED,
            self.CURVE2FLOOD_FLDPLNPY,
            self.CURVE2FLOOD_MULTI_POINT_INTERPOLATION,
        }
    
    def is_curve2flood_fldpln_mapper(self) -> bool:
        return self == self.CURVE2FLOOD_FLDPLNPY