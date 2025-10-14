## Interconvert WGS84 lat/lon, MGRS and UTM coordinate systems on Arduino

This repository hosts a single Arduino program that is a collection of functions from other sources (with links to the originals)
that functions to interconvert GPS coordinate systems WGS84, UTM and MGRS. Latitude is restricted to range from -80 to 84 degrees.

No operational instructions are given here, but the test program should run "out of the box" and demonstrates the basic features. It generates 500 lat/lon pairs uniformly distributed around the Earth, converts each pair to the corresponding MGRS 15 character string, converts that string back to UTM, then converts the UTM back to lat/lon.

Printed output is in .CSV format for spreadsheet evaluation. In my tests with an Arduino Uno R3 (32 bit single precision floats)
the input and output lat/lon pairs differ by at most about 1 meter, so the accuracy is acceptable for normal outdoor activites. I have not yet tested the setup on an Arduino that supports 64 bit double precision. Will update this after doing that.

To check accuracy and verify correctness of MGRS symbols, I use and recommend this all in one online converter

https://rcn.montana.edu/Resources/Converter.aspx

The .h files should be placed in the same sketch folder as the Arduino program.
