#include <GL/glut.h>
#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <math.h>
#include <string.h>




/*pantalla ecran variables*/
int alto = 0,ancho =0;
const char* wndname = "CAMARA";
void salir();
void salirCV();

/*opengl variables*/
const float CAJA_LARGO = 8.0f;
const float CAJA_ANCHO = 2.0f;
const float CAJA_ALTO = 5.0f;
/*openCV variables*/
CvMemStorage* storage = 0;
CvCapture* capture = 0;
IplImage* img = 0;
IplImage* img0 = 0;
int thresh = 100;
CvMat *H_image_vers_modele, *H_modele_vers_image;
float H_GL[16];
float modele[8] = {-CAJA_LARGO/2, -CAJA_ANCHO/2, CAJA_LARGO/2, -CAJA_ANCHO/2, CAJA_LARGO/2, CAJA_ANCHO/2, -CAJA_LARGO/2, CAJA_ANCHO/2};
float amers[8];
