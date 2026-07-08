

class NencartaBaseException(Exception):
    """Base exception for nencarta."""
    pass

class NoStreamsFoundException(NencartaBaseException):
    """Raised when no streams are found in the input geometry."""
    pass
