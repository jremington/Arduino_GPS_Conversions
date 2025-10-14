// Arduino version converts lat/lon to mgrs to utm to lat/lon
// latitude range is -80 to 84
// random lat/lon test OK!
// checked values against this on line WGS84/UTM/MGRS converter: https://rcn.montana.edu/Resources/Converter.aspx

#include "utm.h"
#include "mgrs_to_utm.h"

void setup() {
  Serial.begin(115200);
  while (!Serial) delay(1);
  Serial.println("\n\nTest Suite MGRS/UTM/WGS84\n");

  double lat = 0.0, lon = 0.0;
  long Zone = 0;
  char Hemisphere = 0;
  long Easting = 0;
  long Northing = 0;
  long retval = 0;

  // mgrs test inputs, limited and full precision
  //  char mgrs[] = "49PAN730040"; //=> 10.877072 108.009193 => 49PAN7299803998
  //  char mgrs[] = "10TEP0000071873"; //>= 44.0, -123.0

  char mgrs[20] = {0};

  // generate random lat/lon pairs in range

  for (int i = 0; i < 500; i++) {
    lat = random(-800, 841) / 10.0;
    lon = random(-1800, 1801) / 10.0;

    // Convert WGS84 lat/lon to MGRS

    LLtoMGRS(lat, lon, mgrs);
    Serial.print(lat, 2);  Serial.print(", ");
    Serial.print(lon, 2);  Serial.print(", ");
    Serial.print(mgrs); Serial.print(", ");

    // Convert MGRS to UTM

    retval = Convert_MGRS_To_UTM(mgrs, &Zone, &Hemisphere, &Easting, &Northing);
    if (retval == MGRS_NO_ERROR) {
      Serial.print(Zone); Serial.print(", ");
      Serial.print(Hemisphere); Serial.print(", ");
      Serial.print(Easting); Serial.print(", ");
      Serial.print(Northing); Serial.print(", ");

      //  Convert UTM back to lat/lon

      utmToLatLon (Easting, Northing, Zone, Hemisphere == 'N', lat, lon);
      Serial.print(lat, 6); Serial.print(", ");
      Serial.println(lon, 6);
    }
    else {
      Serial.print("Err #");
      Serial.println(retval);
    }
  }
}


void loop() {}


// ===================== ll to mgrs ==========================
// This from
// https://forums.adafruit.com/viewtopic.php?t=220454
// https://en.wikipedia.org/wiki/Military_Grid_Reference_System

double deg2rad(double d) {
  return d * (M_PI / 180.0);
}
double rad2deg(double r) {
  return r * (180.0 / M_PI);
}

//  valid latitude range -80 to 84 degrees

char UTMLetterDesignator(double Lat) {
  if (Lat < -80) return 'Z';  //error return
  if (Lat >= 72 && Lat <= 84) return 'X';  //special case not in lookup SJR
  const char letters[] = "CDEFGHJKLMNPQRSTUVWX";
  int i = (Lat + 80) / 8;
  return letters[i];
}
// compute digraph of two characters, return as C-string

void MGRSZoneDesignator(double UTMEasting, double UTMNorthing, int zoneNumber, char digraph[]) {
  // Calculate which 100km grid square we're in
  int col = int32_t ((UTMEasting - 100000) / 100000);  // Subtract 100km offset
  int row = int32_t(UTMNorthing / 100000) % 20;

  // MGRS Easting Letters - the pattern shifts every 3 zones
  // Set 0 (zones 1,4,7,10,13,16,19...): ABCDEFGH
  // Set 1 (zones 2,5,8,11,14,17,20...): JKLMNPQR
  // Set 2 (zones 3,6,9,12,15,18,21...): STUVWXYZ

  char eastingLetters[9] = {0};
  int zoneSet = (zoneNumber - 1) % 3;
  if (zoneSet == 0) strcpy(eastingLetters, "ABCDEFGH");
  if (zoneSet == 1) strcpy(eastingLetters, "JKLMNPQR");
  if (zoneSet == 2) strcpy(eastingLetters, "STUVWXYZ");

  // Northing letters are always the same (20 letters, cycling)
  char northingLetters[] = "ABCDEFGHJKLMNPQRSTUV"; // No I or O
  digraph[0] = eastingLetters[col % 8];
  if ((zoneNumber & 1) == 0) row += 5; //start with F in even zones, A in odd for MGRS-new, SJR
  digraph[1] = northingLetters[row % 20];
  digraph[2] = 0; //terminate
}

// convert lat/lon to MGRS string, return in C-string
void LLtoMGRS(double Lat, double Long, char MGRSstring[]) {
  // WGS84 ellipsoid parameters
  const double a = 6378137.0;
  const double f = 1.0 / 298.257223563;
  const double eccSquared = 2 * f - f * f;
  const double k0 = 0.9996;

  double LongTemp = (Long + 180) - int32_t((Long + 180) / 360) * 360 - 180;
  double LatRad = deg2rad(Lat);
  double LongRad = deg2rad(LongTemp);
  uint8_t ZoneNumber = uint8_t((LongTemp + 180) / 6) + 1; //1-60

  // Special cases for Norway and Svalbard
  if (Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0)
    ZoneNumber = 32;
  if ( Lat >= 72.0 && Lat < 84.0 ) {
    if  ( LongTemp >= 0.0  && LongTemp <  9.0 ) ZoneNumber = 31;
    else if ( LongTemp >= 9.0  && LongTemp < 21.0 ) ZoneNumber = 33;
    else if (LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
    else if (LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
  }
  double LongOrigin = (ZoneNumber - 1) * 6 - 180 + 3;
  double LongOriginRad = deg2rad(LongOrigin);

  // More accurate UTM conversion
  double eccPrimeSquared = eccSquared / (1 - eccSquared);
  double N = a / sqrt(1 - eccSquared * sin(LatRad) * sin(LatRad));
  double T = tan(LatRad) * tan(LatRad);
  double C = eccPrimeSquared * cos(LatRad) * cos(LatRad);
  double A = cos(LatRad) * (LongRad - LongOriginRad);

  double M = a * ((1 - eccSquared / 4 - 3 * pow(eccSquared, 2) / 64 - 5 * pow(eccSquared, 3) / 256) * LatRad
                  - (3 * eccSquared / 8 + 3 * pow(eccSquared, 2) / 32 + 45 * pow(eccSquared, 3) / 1024) * sin(2 * LatRad)
                  + (15 * pow(eccSquared, 2) / 256 + 45 * pow(eccSquared, 3) / 1024) * sin(4 * LatRad)
                  - (35 * pow(eccSquared, 3) / 3072) * sin(6 * LatRad));

  double UTMEasting = k0 * N * (A + (1 - T + C) * pow(A, 3) / 6 +
                                (5 - 18 * T + T * T + 72 * C - 58 * eccPrimeSquared) * pow(A, 5) / 120) + 500000.0;
  double UTMNorthing = k0 * (M + N * tan(LatRad) * (A * A / 2 +
                             (5 - T + 9 * C + 4 * C * C) * pow(A, 4) / 24 +
                             (61 - 58 * T + T * T + 600 * C - 330 * eccPrimeSquared) * pow(A, 6) / 720));

  if (Lat < 0) UTMNorthing += 10000000.0;
  // More accurate grid square calculation
  char digraph[3] = {0};
  MGRSZoneDesignator(UTMEasting, UTMNorthing, ZoneNumber, digraph);

  int32_t easting = int32_t(UTMEasting) % 100000;
  int32_t northing = int32_t(UTMNorthing) % 100000;

  sprintf(MGRSstring, "%02d%c%s%05ld%05ld", ZoneNumber, UTMLetterDesignator(Lat), digraph, easting, northing);
  return;
}

// mgrs_to_utm code taken from
// https://www.stellman-greene.com/mgrs_to_utm/
// Ellipsoid parameters, default to WGS 84

typedef struct Latitude_Band_Value
{
  long letter;            /* letter representing latitude band  */
  double min_northing;    /* minimum northing for latitude band */
  double north;           /* upper latitude for latitude band   */
  double south;           /* lower latitude for latitude band   */
} Latitude_Band;

const Latitude_Band Latitude_Band_Table[20] =
{ {LETTER_C, 1100000.0, -72.0, -80.5},
  {LETTER_D, 2000000.0, -64.0, -72.0},
  {LETTER_E, 2800000.0, -56.0, -64.0},
  {LETTER_F, 3700000.0, -48.0, -56.0},
  {LETTER_G, 4600000.0, -40.0, -48.0},
  {LETTER_H, 5500000.0, -32.0, -40.0},
  {LETTER_J, 6400000.0, -24.0, -32.0},
  {LETTER_K, 7300000.0, -16.0, -24.0},
  {LETTER_L, 8200000.0, -8.0, -16.0},
  {LETTER_M, 9100000.0, 0.0, -8.0},
  {LETTER_N, 0.0, 8.0, 0.0},
  {LETTER_P, 800000.0, 16.0, 8.0},
  {LETTER_Q, 1700000.0, 24.0, 16.0},
  {LETTER_R, 2600000.0, 32.0, 24.0},
  {LETTER_S, 3500000.0, 40.0, 32.0},
  {LETTER_T, 4400000.0, 48.0, 40.0},
  {LETTER_U, 5300000.0, 56.0, 48.0},
  {LETTER_V, 6200000.0, 64.0, 56.0},
  {LETTER_W, 7000000.0, 72.0, 64.0},
  {LETTER_X, 7900000.0, 84.5, 72.0}
};

long Convert_MGRS_To_UTM (char   *MGRS,
                          long   * Zone,
                          char   *Hemisphere,
                          long * Easting,
                          long * Northing)
/*
   The function Convert_MGRS_To_UTM converts an MGRS coordinate string
   to UTM projection (zone, hemisphere, easting and northing) coordinates
   according to the current ellipsoid parameters.  If any errors occur,
   the error code(s) are returned by the function, otherwise UTM_NO_ERROR
   is returned.

      MGRS       : MGRS coordinate string           (input)
      Zone       : UTM zone                         (output)
      Hemisphere : North or South hemisphere        (output)
      Easting    : Easting (X) in meters            (output)
      Northing   : Northing (Y) in meters           (output)
*/

{ /* Convert_MGRS_To_UTM */
  double scaled_min_northing;
  double min_northing;
  long ltr2_low_value;
  long ltr2_high_value;
  double false_northing;
  double grid_easting;        /* Easting for 100,000 meter grid square      */
  double grid_northing;       /* Northing for 100,000 meter grid square     */
  uint8_t letters[3] = {0};
  long in_precision;
  long error_code = MGRS_NO_ERROR;

  error_code = Break_MGRS_String (MGRS, Zone, letters, Easting, Northing, &in_precision);
  if (!*Zone)
    error_code |= MGRS_STRING_ERROR;
  else
  {
    if (!error_code)
    {
      if ((letters[0] == LETTER_X) && ((*Zone == 32) || (*Zone == 34) || (*Zone == 36)))
        error_code |= MGRS_STRING_ERROR;
      else
      {
        if (letters[0] < LETTER_N)
          *Hemisphere = 'S';
        else
          *Hemisphere = 'N';

        Get_Grid_Values(*Zone, &ltr2_low_value, &ltr2_high_value, &false_northing);

        /* Check that the second letter of the MGRS string is within
           the range of valid second letter values
           Also check that the third letter is valid */
        if ((letters[1] < ltr2_low_value) || (letters[1] > ltr2_high_value) || (letters[2] > LETTER_V))
          error_code |= MGRS_STRING_ERROR;

        if (!error_code)
        {
          grid_northing = (double)(letters[2]) * ONEHT + false_northing;
          grid_easting = (double)((letters[1]) - ltr2_low_value + 1) * ONEHT;
          if ((ltr2_low_value == LETTER_J) && (letters[1] > LETTER_O))
            grid_easting = grid_easting - ONEHT;

          if (letters[2] > LETTER_O)
            grid_northing = grid_northing - ONEHT;

          if (letters[2] > LETTER_I)
            grid_northing = grid_northing - ONEHT;

          if (grid_northing >= TWOMIL)
            grid_northing = grid_northing - TWOMIL;

          error_code = Get_Latitude_Band_Min_Northing(letters[0], &min_northing);
          if (!error_code)
          {
            scaled_min_northing = min_northing;
            while (scaled_min_northing >= TWOMIL)
            {
              scaled_min_northing = scaled_min_northing - TWOMIL;
            }

            grid_northing = grid_northing - scaled_min_northing;
            if (grid_northing < 0.0)
              grid_northing = grid_northing + TWOMIL;

            grid_northing = min_northing + grid_northing;

            *Easting = (long) grid_easting + *Easting;
            *Northing = (long) grid_northing + *Northing;
          }
        }
      }
    }
  }
  return (error_code);
} /* END Convert_MGRS_To_UTM */


long Break_MGRS_String (char* MGRS,
                        long * Zone,
                        uint8_t Letters[3],
                        long * Easting,
                        long * Northing,
                        long * Precision)
/*
   The function Break_MGRS_String breaks down an MGRS
   coordinate string into its component parts.

     MGRS           : MGRS coordinate string          (input)
     Zone           : UTM Zone                        (output)
     Letters        : MGRS coordinate string letters  (output)
     Easting        : Easting value                   (output)
     Northing       : Northing value                  (output)
     Precision      : Precision level of MGRS string  (output)
*/
{ /* Break_MGRS_String */
  long num_digits;
  long num_letters;
  long i = 0;
  long j = 0;
  long error_code = MGRS_NO_ERROR;

  while (MGRS[i] == ' ')
    i++;  /* skip any leading blanks */
  j = i;
  while (isdigit(MGRS[i]))
    i++;
  num_digits = i - j;
  if (num_digits <= 2)
    if (num_digits > 0)
    {
      char zone_string[3] = {0};
      /* get zone */
      strncpy (zone_string, MGRS + j, 2);
      zone_string[2] = 0;
      *Zone = atol(zone_string);
      if ((*Zone < 1) || (*Zone > 60))
        error_code |= MGRS_STRING_ERROR;
    }
    else
      *Zone = 0;
  else
    error_code |= MGRS_STRING_ERROR;
  j = i;

  while (isalpha(MGRS[i]))
    i++;
  num_letters = i - j;
  if (num_letters == 3)
  {
    /* get letters */
    Letters[0] = (toupper(MGRS[j]) - 'A');
    if ((Letters[0] == LETTER_I) || (Letters[0] == LETTER_O))
      error_code |= MGRS_STRING_ERROR;
    Letters[1] = (toupper(MGRS[j + 1]) - 'A');
    if ((Letters[1] == LETTER_I) || (Letters[1] == LETTER_O))
      error_code |= MGRS_STRING_ERROR;
    Letters[2] = (toupper(MGRS[j + 2]) - 'A');
    if ((Letters[2] == LETTER_I) || (Letters[2] == LETTER_O))
      error_code |= MGRS_STRING_ERROR;
  }
  else
    error_code |= MGRS_STRING_ERROR;
  j = i;
  while (isdigit(MGRS[i]))
    i++;
  num_digits = i - j;
  if ((num_digits <= 10) && (num_digits % 2 == 0))
  {
    long n;
    char east_string[6];
    char north_string[6];
    long east;
    long north;
    double multiplier;
    /* get easting & northing */
    n = num_digits / 2;
    *Precision = n;
    if (n > 0)
    {
      strncpy (east_string, MGRS + j, n);
      east_string[n] = 0;
      east = atol(east_string);
      strncpy (north_string, MGRS + j + n, n);
      north_string[n] = 0;
      north = atol(north_string);
      multiplier = pow (10.0, 5 - n);
      *Easting = (long) (east * multiplier);
      *Northing = (long) (north * multiplier);
    }
    else
    {
      *Easting = 0;
      *Northing = 0;
    }
  }
  else
    error_code |= MGRS_STRING_ERROR;

  return (error_code);
} /* Break_MGRS_String */



long Get_Latitude_Band_Min_Northing(long letter, double * min_northing)
/*
   The function Get_Latitude_Band_Min_Northing receives a latitude band letter
   and uses the Latitude_Band_Table to determine the minimum northing for that
   latitude band letter.

     letter        : Latitude band letter             (input)
     min_northing  : Minimum northing for that letter(output)
*/
{ /* Get_Latitude_Band_Min_Northing */
  long error_code = MGRS_NO_ERROR;

  if ((letter >= LETTER_C) && (letter <= LETTER_H))
    *min_northing = Latitude_Band_Table[letter - 2].min_northing;
  else if ((letter >= LETTER_J) && (letter <= LETTER_N))
    *min_northing = Latitude_Band_Table[letter - 3].min_northing;
  else if ((letter >= LETTER_P) && (letter <= LETTER_X))
    *min_northing = Latitude_Band_Table[letter - 4].min_northing;
  else
    error_code |= MGRS_STRING_ERROR;

  return error_code;
} /* END Get_Latitude_Band_Min_Northing */


void Get_Grid_Values (long zone,
                      long * ltr2_low_value,
                      long * ltr2_high_value,
                      double * false_northing)
/*
   The function Get_Grid_Values sets the letter range used for
   the 2nd letter in the MGRS coordinate string, based on the set
   number of the utm zone. It also sets the false northing using a
   value of A for the second letter of the grid square, based on
   the grid pattern and set number of the utm zone.

      zone            : Zone number             (input)
      ltr2_low_value  : 2nd letter low number   (output)
      ltr2_high_value : 2nd letter high number  (output)
      false_northing  : False northing          (output)
*/
{ /* BEGIN Get_Grid_Values */
  long set_number;    /* Set number (1-6) based on UTM zone number */
  set_number = zone % 6;

  if (!set_number)
    set_number = 6;

  if ((set_number == 1) || (set_number == 4))
  {
    *ltr2_low_value = LETTER_A;
    *ltr2_high_value = LETTER_H;
  }
  else if ((set_number == 2) || (set_number == 5))
  {
    *ltr2_low_value = LETTER_J;
    *ltr2_high_value = LETTER_R;
  }
  else if ((set_number == 3) || (set_number == 6))
  {
    *ltr2_low_value = LETTER_S;
    *ltr2_high_value = LETTER_Z;
  }

  /* False northing at A for second letter of grid square */

  if ((set_number % 2) ==  0)
    *false_northing = 1500000.0;
  else
    *false_northing = 0.0;

} /* END OF Get_Grid_Values */

// UTM.c  convert from WGS84 lat/long to UTM
// (the code for the reverse operation from this source fails on some inputs)
// https://alephnull.net/software/gis/UTM_WGS84_C_plus_plus.shtml
// Original Javascript by Chuck Taylor
// Port to C++ by Alex Hajnal
//
// *** THIS CODE USES 32-BIT FLOATS BY DEFAULT ***
// *** For 64-bit double-precision edit UTM.h: undefine FLOAT_32 and define FLOAT_64
//
// This is a simple port of the code on the Geographic/UTM Coordinate Converter (1) page from Javascript to C++.
// Using this you can easily convert between UTM and WGS84 (latitude and longitude).
// Accuracy seems to be around 50cm (I suspect rounding errors are limiting precision).
// This code is provided as-is and has been minimally tested; enjoy but use at your own risk!
// The license for UTM.cpp and UTM.h is the same as the original Javascript:
// "The C++ source code in UTM.cpp and UTM.h may be copied and reused without restriction."


// DegToRad
// Converts degrees to radians.
FLOAT DegToRad(FLOAT deg) {
  return (M_PI * deg / 180.0);
}


// RadToDeg
// Converts radians to degrees.
FLOAT RadToDeg(FLOAT rad) {
  return (rad / M_PI * 180.0);
}

// ArcLengthOfMeridian
// Computes the ellipsoidal distance from the equator to a point at a
// given latitude.
//
// Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
// GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
//
// Inputs:
//     phi - Latitude of the point, in radians.
//
// Globals:
//     sm_a - Ellipsoid model major axis.
//     sm_b - Ellipsoid model minor axis.
//
// Returns:
//     The ellipsoidal distance of the point from the equator, in meters.
FLOAT ArcLengthOfMeridian (FLOAT phi) {
  FLOAT alpha, beta, gamma, delta, epsilon, n;
  FLOAT result;

  /* Precalculate n */
  n = (sm_a - sm_b) / (sm_a + sm_b);

  /* Precalculate alpha */
  alpha = ((sm_a + sm_b) / 2.0)
          * (1.0 + (POW(n, 2.0) / 4.0) + (POW(n, 4.0) / 64.0));

  /* Precalculate beta */
  beta = (-3.0 * n / 2.0) + (9.0 * POW(n, 3.0) / 16.0)
         + (-3.0 * POW(n, 5.0) / 32.0);

  /* Precalculate gamma */
  gamma = (15.0 * POW(n, 2.0) / 16.0)
          + (-15.0 * POW(n, 4.0) / 32.0);

  /* Precalculate delta */
  delta = (-35.0 * POW(n, 3.0) / 48.0)
          + (105.0 * POW(n, 5.0) / 256.0);

  /* Precalculate epsilon */
  epsilon = (315.0 * POW(n, 4.0) / 512.0);

  /* Now calculate the sum of the series and return */
  result = alpha
           * (phi + (beta * SIN(2.0 * phi))
              + (gamma * SIN(4.0 * phi))
              + (delta * SIN(6.0 * phi))
              + (epsilon * SIN(8.0 * phi)));

  return result;
}

// UTMCentralMeridian
// Determines the central meridian for the given UTM zone.
//
// Inputs:
//     zone - An integer value designating the UTM zone, range [1,60].
//
// Returns:
//   The central meridian for the given UTM zone, in radians
//   Range of the central meridian is the radian equivalent of [-177,+177].
FLOAT UTMCentralMeridian(int zone) {
  FLOAT cmeridian;
  cmeridian = DegToRad(-183.0 + ((FLOAT)zone * 6.0));

  return cmeridian;
}

// FootpointLatitude
//
// Computes the footpoint latitude for use in converting transverse
// Mercator coordinates to ellipsoidal coordinates.
//
// Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
//   GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
//
// Inputs:
//   y - The UTM northing coordinate, in meters.
//
// Returns:
//   The footpoint latitude, in radians.
FLOAT FootpointLatitude(FLOAT y) {
  FLOAT y_, alpha_, beta_, gamma_, delta_, epsilon_, n;
  FLOAT result;

  /* Precalculate n (Eq. 10.18) */
  n = (sm_a - sm_b) / (sm_a + sm_b);

  /* Precalculate alpha_ (Eq. 10.22) */
  /* (Same as alpha in Eq. 10.17) */
  alpha_ = ((sm_a + sm_b) / 2.0)
           * (1 + (POW(n, 2.0) / 4) + (POW(n, 4.0) / 64));

  /* Precalculate y_ (Eq. 10.23) */
  y_ = y / alpha_;

  /* Precalculate beta_ (Eq. 10.22) */
  beta_ = (3.0 * n / 2.0) + (-27.0 * POW(n, 3.0) / 32.0)
          + (269.0 * POW(n, 5.0) / 512.0);

  /* Precalculate gamma_ (Eq. 10.22) */
  gamma_ = (21.0 * POW(n, 2.0) / 16.0)
           + (-55.0 * POW(n, 4.0) / 32.0);

  /* Precalculate delta_ (Eq. 10.22) */
  delta_ = (151.0 * POW(n, 3.0) / 96.0)
           + (-417.0 * POW(n, 5.0) / 128.0);

  /* Precalculate epsilon_ (Eq. 10.22) */
  epsilon_ = (1097.0 * POW(n, 4.0) / 512.0);

  /* Now calculate the sum of the series (Eq. 10.21) */
  result = y_ + (beta_ * SIN(2.0 * y_))
           + (gamma_ * SIN(4.0 * y_))
           + (delta_ * SIN(6.0 * y_))
           + (epsilon_ * SIN(8.0 * y_));

  return result;
}

// MapLatLonToXY
// Converts a latitude/longitude pair to x and y coordinates in the
// Transverse Mercator projection.  Note that Transverse Mercator is not
// the same as UTM; a scale factor is required to convert between them.
//
// Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
// GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
//
// Inputs:
//    phi - Latitude of the point, in radians.
//    lambda - Longitude of the point, in radians.
//    lambda0 - Longitude of the central meridian to be used, in radians.
//
// Outputs:
//    x - The x coordinate of the computed point, flaat or double
//    y - The y coordinate of the computed point, float or double
//
// Returns:
//    The function does not return a value.
void MapLatLonToXY (FLOAT phi, FLOAT lambda, FLOAT lambda0, FLOAT & x, FLOAT & y) {
  FLOAT N, nu2, ep2, t, t2, l;
  FLOAT l3coef, l4coef, l5coef, l6coef, l7coef, l8coef;

  /* Precalculate ep2 */
  ep2 = (sm_a * sm_a - sm_b * sm_b) / (sm_b * sm_b);

  /* Precalculate nu2 */
  nu2 = ep2 * COS(phi) * COS(phi);

  /* Precalculate N */
  N = sm_a * sm_a / (sm_b * SQRT(1 + nu2));

  /* Precalculate t */
  t = TAN(phi);
  t2 = t * t;

  /* Precalculate l */
  l = lambda - lambda0;

  /* Precalculate coefficients for l**n in the equations below
     so a normal human being can read the expressions for easting
     and northing
     -- l**1 and l**2 have coefficients of 1.0 */
  l3coef = 1.0 - t2 + nu2;

  l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);

  l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2
           - 58.0 * t2 * nu2;

  l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2
           - 330.0 * t2 * nu2;

  l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);

  l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);

  /* Calculate easting (x) */
  x = N * COS(phi) * l
      + (N / 6.0 * POW(COS(phi), 3.0) * l3coef * POW(l, 3.0))
      + (N / 120.0 * POW(COS(phi), 5.0) * l5coef * POW(l, 5.0))
      + (N / 5040.0 * POW(COS(phi), 7.0) * l7coef * POW(l, 7.0));

  /* Calculate northing (y) */
  y = ArcLengthOfMeridian (phi)
      + (t / 2.0 * N * POW(COS(phi), 2.0) * POW(l, 2.0))
      + (t / 24.0 * N * POW(COS(phi), 4.0) * l4coef * POW(l, 4.0))
      + (t / 720.0 * N * POW(COS(phi), 6.0) * l6coef * POW(l, 6.0))
      + (t / 40320.0 * N * POW(COS(phi), 8.0) * l8coef * POW(l, 8.0));

  return;
}

// LatLonToUTMXY
// Converts a latitude/longitude pair to x and y coordinates in the
// Universal Transverse Mercator projection.
//
// Inputs:
//   lat - Latitude of the point, in degrees
//   lon - Longitude of the point, in degrees.
//   zone - UTM zone to be used for calculating values for x and y.
//          If zone is less than 1 or greater than 60, the routine
//          will determine the appropriate zone from the value of lon.
//
// Outputs:
//   x - The x coordinate (easting) of the computed point. (in meters)
//   y - The y coordinate (northing) of the computed point. (in meters)
//
// Returns:
//   The UTM zone used for calculating the values of x and y.

int LatLonToUTMXY (FLOAT lat, FLOAT lon, int zone, long & x, long & y) {
  if ( (zone < 1) || (zone > 60) )
    zone = FLOOR((lon + 180.0) / 6) + 1;
  FLOAT xx = x;
  FLOAT yy = y;
  MapLatLonToXY (DegToRad(lat), DegToRad(lon), UTMCentralMeridian(zone), xx, yy);

  /* Adjust easting and northing for UTM system. */
  xx = xx * UTMScaleFactor + 500000.0;
  yy = yy * UTMScaleFactor;
  if (yy < 0.0)
    yy = yy + 10000000.0;
  x = (long) xx;
  y = (long) yy;
  return zone;
}

// Convert UTM to WGS84 lat/lon
// from  https://github.com/tecfinite/utm_convert/blob/main/main.cpp
// added error return, long ints for easting and northing instead of double

int utmToLatLon(long easting, long northing, int zone, bool isNorthernHemisphere, double & lat, double & lon) {
  double a = 6378137.0;
  double f = 1.0 / 298.257223563;
  double k0 = 0.9996;
  double e = sqrt(f * (2 - f));
  double e1sq = e * e / (1 - e * e);
  if (zone < 1 || zone > 60) {
    lat = 0.0;
    lon = 0.0;
    return 1; //zone error
  }
  double x = easting - 500000.0;
  double y = northing;
  if (!isNorthernHemisphere) {
    y -= 10000000.0;
  }

  double M = y / k0;
  double mu = M / (a * (1 - e * e / 4 - 3 * pow(e, 4) / 64 - 5 * pow(e, 6) / 256));

  double e1 = (1 - sqrt(1 - e * e)) / (1 + sqrt(1 - e * e));

  double phi1Rad = mu
                   + (3 * e1 / 2 - 27 * pow(e1, 3) / 32) * sin(2 * mu)
                   + (21 * pow(e1, 2) / 16 - 55 * pow(e1, 4) / 32) * sin(4 * mu)
                   + (151 * pow(e1, 3) / 96) * sin(6 * mu);

  double N1 = a / sqrt(1 - pow(e * sin(phi1Rad), 2));
  double T1 = pow(tan(phi1Rad), 2);
  double C1 = e1sq * pow(cos(phi1Rad), 2);
  double R1 = a * (1 - e * e) / pow(1 - pow(e * sin(phi1Rad), 2), 1.5);
  double D = x / (N1 * k0);

  double latRad = phi1Rad - (N1 * tan(phi1Rad) / R1) *
                  (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * e1sq) * pow(D, 4) / 24
                   + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * e1sq - 3 * C1 * C1) * pow(D, 6) / 720);
  lat = latRad * 180.0 / PI;

  double lonOrigin = (zone - 1) * 6.0 - 180.0 + 3.0;
  double lonRad = (D - (1 + 2 * T1 + C1) * pow(D, 3) / 6
                   + (5 - 2 * C1 + 28 * T1 - 3 * pow(C1, 2) + 8 * e1sq + 24 * pow(T1, 2)) * pow(D, 5) / 120)
                  / cos(phi1Rad);
  lon = lonOrigin + lonRad * 180.0 / PI;
  return 0;
}
