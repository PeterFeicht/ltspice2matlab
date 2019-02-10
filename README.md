# LTspice2Matlab

Load .raw simulation files created with LTspice into MATLAB.

**LTspice2Matlab** imports an LTspice IV or LTspice XVII .raw waveform file
containing data from Transient Analysis (.tran), AC Analysis (.ac), DC Sweep
(.dc), Operating Point (.op), Transfer Function (.tf), FFT (.four), or Noise
(.noise) simulation, and converts voltages and currents vs. time (or frequency)
into a MATLAB data structure. It also supports import of stepped simulations.
This function can read compressed binary, uncompressed binary, and ASCII file
formats. In the case of compressed binary, the data is automatically
uncompressed using fast quadratic point insertion. This function handles very
large binary simulation files efficiently, and has an option to load only a
subset of a file's waveforms to reduce memory consumption. Type
`>> help LTspice2Matlab` for details.

Use **LTspice2Matlab** to import LTspice waveforms for additional analysis in
MATLAB, or for comparison with measured data. **LTspice2Matlab** has been
tested with LTspice IV and LTspice XVII, and supports MATLAB versions 2016b
and later (support for earlier MATLAB version may be restored by replacing a
few functions with equivalent, but less efficient, alternatives).

The original version of this script created by Paul Wagner can be found on
[MATLAB Central][original].

[original]: https://www.mathworks.com/matlabcentral/fileexchange/23394-fast-import-of-compressed-binary-raw-files-created-with-ltspice-circuit-simulator
