/* mgrs_to_utm - convert MGRS coordinates to UTM
 *
 * http://svn.stellman-greene.com/mgrs_to_utm/trunk/
 *
 * to install and run: 
 *   - download mgrs_to_utm.c and mgrs_to_utm.h
 *   - build: gcc -o mgrs_to_utm mgrs_to_utm.c
 *   - run mgrs_to_utm
 *
 * For more information on MGRS and UTM coordinates:
 *   Military Grid Reference System 
 *      http://en.wikipedia.org/wiki/Military_grid_reference_system
 *   Universal Transverse Mercator coordinate system
 *      http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
 *
 * This project was adapted from Geotrans 2.4.1
 *    http://earth-info.nga.mil/GandG/geotrans/
 *
 */

#ifndef MGRS_TO_UTM_H
#define MGRS_TO_UTM_H

#define MGRS_NO_ERROR                0x0000
#define MGRS_LAT_ERROR               0x0001
#define MGRS_LON_ERROR               0x0002
#define MGRS_STRING_ERROR            0x0004
#define MGRS_PRECISION_ERROR         0x0008
#define MGRS_A_ERROR                 0x0010
#define MGRS_INV_F_ERROR             0x0020
#define MGRS_EASTING_ERROR           0x0040
#define MGRS_NORTHING_ERROR          0x0080
#define MGRS_ZONE_ERROR              0x0100
#define MGRS_HEMISPHERE_ERROR        0x0200


#define LETTER_A               0   /* ARRAY INDEX FOR LETTER A               */
#define LETTER_B               1   /* ARRAY INDEX FOR LETTER B               */
#define LETTER_C               2   /* ARRAY INDEX FOR LETTER C               */
#define LETTER_D               3   /* ARRAY INDEX FOR LETTER D               */
#define LETTER_E               4   /* ARRAY INDEX FOR LETTER E               */
#define LETTER_F               5   /* ARRAY INDEX FOR LETTER E               */
#define LETTER_G               6   /* ARRAY INDEX FOR LETTER H               */
#define LETTER_H               7   /* ARRAY INDEX FOR LETTER H               */
#define LETTER_I               8   /* ARRAY INDEX FOR LETTER I               */
#define LETTER_J               9   /* ARRAY INDEX FOR LETTER J               */
#define LETTER_K              10   /* ARRAY INDEX FOR LETTER J               */
#define LETTER_L              11   /* ARRAY INDEX FOR LETTER L               */
#define LETTER_M              12   /* ARRAY INDEX FOR LETTER M               */
#define LETTER_N              13   /* ARRAY INDEX FOR LETTER N               */
#define LETTER_O              14   /* ARRAY INDEX FOR LETTER O               */
#define LETTER_P              15   /* ARRAY INDEX FOR LETTER P               */
#define LETTER_Q              16   /* ARRAY INDEX FOR LETTER Q               */
#define LETTER_R              17   /* ARRAY INDEX FOR LETTER R               */
#define LETTER_S              18   /* ARRAY INDEX FOR LETTER S               */
#define LETTER_T              19   /* ARRAY INDEX FOR LETTER S               */
#define LETTER_U              20   /* ARRAY INDEX FOR LETTER U               */
#define LETTER_V              21   /* ARRAY INDEX FOR LETTER V               */
#define LETTER_W              22   /* ARRAY INDEX FOR LETTER W               */
#define LETTER_X              23   /* ARRAY INDEX FOR LETTER X               */
#define LETTER_Y              24   /* ARRAY INDEX FOR LETTER Y               */
#define LETTER_Z              25   /* ARRAY INDEX FOR LETTER Z               */
#define MGRS_LETTERS            3  /* NUMBER OF LETTERS IN MGRS              */
#define ONEHT          100000.e0    /* ONE HUNDRED THOUSAND                  */
#define TWOMIL        2000000.e0    /* TWO MILLION                           */
#define TRUE                      1  /* CONSTANT VALUE FOR TRUE VALUE  */
#define FALSE                     0  /* CONSTANT VALUE FOR FALSE VALUE */
//#define PI    3.14159265358979323e0  /* PI                             */
#define PI_OVER_2  (PI / 2.0e0)
/*

long Convert_MGRS_To_UTM (char   *MGRS,
			  long   *Zone,
			  char   *Hemisphere,
			  double *Easting,
			  double *Northing); 


long Break_MGRS_String (char* MGRS,
                        long* Zone,
                        long Letters[3],
                        double* Easting,
                        double* Northing,
                        long* Precision);


long Get_Latitude_Band_Min_Northing(long letter, double* min_northing);

void Get_Grid_Values (long zone, 
                      long* ltr2_low_value, 
                      long* ltr2_high_value, 
                      double *false_northing);
*/
#endif
