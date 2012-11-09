# LTspice2Matlab

Loads .raw simulation files created with LTspice.

**LTspice2Matlab** imports an LTspice IV .RAW waveform file containing data
from a Transient Analysis (.tran) or AC Analysis (.ac) simulation, and converts
voltages and currents vs. time (or frequency) into a MATLAB data structure.
This function can read compressed binary, uncompressed binary, and ASCII file
formats. In the case of compressed binary, the data is automatically
uncompressed using fast quadratic point insertion. This function handles very
large binary simulation files efficiently, and has an option to load only a
subset of a file's waveforms to reduce memory consumption. Type
`>> help LTspice2Matlab` for details.

Use **LTspice2Matlab** to import LTspice waveforms for additional analysis in
MATLAB, or for comparison with measured data. **LTspice2Matlab** has been
tested with LTspice IV version 4.01p, and MATLAB versions 6.1 and 7.5, and
regression testing has been used to expose the function to a wide range of
LTspice settings.
