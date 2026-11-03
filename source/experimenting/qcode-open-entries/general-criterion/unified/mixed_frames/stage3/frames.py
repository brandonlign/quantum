"""Stage-3 frames; stage-1 code is imported read-only."""
from functools import lru_cache
from mixed_frames.frames import adapted_frame, mixed_frame as original_frame

@lru_cache(None)
def mixed_frame(name):
    return adapted_frame(8, [0xff00, 0x00ff], name) if name == '888' else original_frame(name)
