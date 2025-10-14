## Interconvert WGS84 lat/lon, MGRS and UTM coordinate systems on Arduino

This repository hosts a single Arduino program consisting of a collection of functions from several sources (with links to the originals)
that operates to interconvert GPS location data, assuming coordinate systems WGS84, UTM and MGRS. Latitude is restricted to range from -80 to 84 degrees.

No operational instructions are given here, but the test program should run "out of the box" and demonstrates the basic features. It generates 500 lat/lon pairs uniformly distributed around the Earth, converts each pair to the corresponding MGRS 15 character string, converts that string back to UTM, then converts the UTM back to lat/lon.

Printed output of the test suite is in .CSV format for spreadsheet evaluation. With an Arduino Uno R3 (32 bit single precision floats)
the input and output lat/lon pairs differ by at most about 1 meter, so the accuracy is acceptable for normal outdoor activites. Also successfully tested on an Adafruit 
Feather M0 (64 bit doubles) with the expected improved accuracy.

To check accuracy and verify correctness of MGRS symbols, I use and recommend this all in one online converter

https://rcn.montana.edu/Resources/Converter.aspx

The .h files should be placed in the same sketch folder as the Arduino program.
