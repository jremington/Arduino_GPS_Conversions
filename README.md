## Interconvert WGS84 lat/lon, MGRS and UTM coordinate systems on Arduino

This repository hosts a single Arduino program consisting of a collection of functions (a library of sorts with links to the original sources),
for interconversion of GPS location data. Coordinate systems WGS84, UTM and MGRS are supported, with latitude restricted to the range from -80 to 84 degrees.

No operational instructions are given here, but the test program should run "out of the box" and demonstrates the basic features. It generates 500 lat/lon pairs uniformly distributed around the Earth, converts each pair to the corresponding MGRS 15 character string, converts that string back to UTM, then converts the UTM back to lat/lon.

Printed output of the test suite is in .CSV format for spreadsheet evaluation. Tested with an Arduino Uno R3 (32 bit single precision floats) and an Adafruit 
Feather SAMD21 M0 (64 bit doubles). 

With 32 bit floats the input and output lat/lon pairs differ by at most about 1 meter, so even in that case the accuracy is comparable to consumer GPS units and acceptable for normal outdoor activities.

To check conversion accuracy and verify correctness of MGRS symbols, I used and recommend this online converter

https://rcn.montana.edu/Resources/Converter.aspx

The .h files should be placed in the same sketch folder as the Arduino program.
