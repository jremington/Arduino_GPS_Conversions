## Interconvert WGS84 lat/lon, MGRS and UTM coordinate systems on Arduino

This repository hosts a single Arduino program that is a collection of functions from other sources (with links to the originals)
that functions interconvert GPS coordinate systems WGS84, UTM and MGRS. Latitude range is restricted to -80 to 84 degrees.

No instructions are given, but the test program shows all the basic features. It generates 500 lat/lon pairs uniformly distributed around the Earth,
converts each pair to the corresponding MGRS 15 character string, converts that string back to UTM, then converts the UTM back to lat/lon.

Printed output is in .CSV format for spreadsheet evaluation. In my tests with an Arduino Uno (32 bit floats)
the input and output lat/lon pairs differ by at most about 1 meter, so the accuracy is acceptable for normal outdoor activites.

To check accuracy and verify correction of MGRS symbols, I use and recommend this all in one online converter

https://rcn.montana.edu/Resources/Converter.aspx
