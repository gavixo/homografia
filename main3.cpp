#include "variables.h"
//#include "gl_func.cpp"
//#include "cv_func.cpp"
void calcul_homographie(CvMat *H_I_vers_M, CvMat *H_M_vers_I, float *H_GL) {
	CvMat src_points = cvMat(4, 2, CV_32F, amers);
	CvMat dst_points = cvMat(4, 2, CV_32F, modele);
	
	cvFindHomography( &src_points, &dst_points, H_I_vers_M );
	cvInv(H_I_vers_M, H_M_vers_I, CV_LU);
	int i, j; float r;
	
	for (j=0; j<H_M_vers_I->cols; j++)
		for (i=0; i<H_M_vers_I->rows; i++) {
			r = cvmGet(H_M_vers_I, j, i)/cvmGet(H_M_vers_I, 2, 2);
			cvmSet(H_M_vers_I, j, i, r);
		}
		
	H_GL[0] = cvmGet(H_M_vers_I, 0, 0); H_GL[4] = cvmGet(H_M_vers_I, 0, 1); H_GL[8] = 0; H_GL[12] = cvmGet(H_M_vers_I, 0, 2);
	H_GL[1] = cvmGet(H_M_vers_I, 1, 0); H_GL[5] = cvmGet(H_M_vers_I, 1, 1); H_GL[9] = 0; H_GL[13] = cvmGet(H_M_vers_I, 1, 2);
	H_GL[2] = 0; H_GL[6] = 0; H_GL[10] = 0; H_GL[14] = 0;
	H_GL[3] = cvmGet(H_M_vers_I, 2, 0); H_GL[7] = cvmGet(H_M_vers_I, 2, 1); H_GL[11] = cvmGet(H_M_vers_I, 2, 2); H_GL[15] = 1;
}
// helper function:
// finds a cosine of angle between vectors
// from pt0->pt1 and from pt0->pt2 
double angle( CvPoint* pt1, CvPoint* pt2, CvPoint* pt0 )
{
    double dx1 = pt1->x - pt0->x;
    double dy1 = pt1->y - pt0->y;
    double dx2 = pt2->x - pt0->x;
    double dy2 = pt2->y - pt0->y;
    return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}
// returns sequence of squares detected on the image.
// the sequence is stored in the specified memory storage
CvSeq* findSquares4( IplImage* img, CvMemStorage* storage )
{
    CvSeq* contours;
    int i, c, l, N = 11;
    CvSize sz = cvSize( img->width & -2, img->height & -2 );
    IplImage* timg = cvCloneImage( img ); // make a copy of input image
    IplImage* gray = cvCreateImage( sz, 8, 1 ); 
    IplImage* pyr = cvCreateImage( cvSize(sz.width/2, sz.height/2), 8, 3 );
    IplImage* tgray;
    CvSeq* result;
    double s, t;
    // create empty sequence that will contain points -
    // 4 points per square (the square's vertices)
    CvSeq* squares = cvCreateSeq( 0, sizeof(CvSeq), sizeof(CvPoint), storage );
    
    // select the maximum ROI in the image
    // with the width and height divisible by 2
    cvSetImageROI( timg, cvRect( 0, 0, sz.width, sz.height ));
    
    // down-scale and upscale the image to filter out the noise
    cvPyrDown( timg, pyr, 7 );
    cvPyrUp( pyr, timg, 7 );
    tgray = cvCreateImage( sz, 8, 1 );
    
    // find squares in every color plane of the image
    for( c = 0; c < 3; c++ )
    {
        // extract the c-th color plane
        cvSetImageCOI( timg, c+1 );
        cvCopy( timg, tgray, 0 );
        
        // try several threshold levels
        for( l = 0; l < N; l++ )
        {
            // hack: use Canny instead of zero threshold level.
            // Canny helps to catch squares with gradient shading   
            if( l == 0 )
            {
                // apply Canny. Take the upper threshold from slider
                // and set the lower to 0 (which forces edges merging) 
                cvCanny( tgray, gray, 10, thresh, 3 );
                //cvShowImage( "d",gray);
                // dilate canny output to remove potential
                // holes between edge segments 
                cvDilate( gray, gray, 0, 1 ); 
                cvShowImage( "d",gray);               
            }
            else
            {
                // apply threshold if l!=0:
                //     tgray(x,y) = gray(x,y) < (l+1)*255/N ? 255 : 0
                cvThreshold( tgray, gray, (l+1)*255/N, 255, CV_THRESH_BINARY );
                 
            }
            
            // find contours and store them all as a list
            int contor = cvFindContours( gray, storage, &contours, sizeof(CvContour),
                CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0,0) );
            
            // test each contour
            while( contours )
            {
                // approximate contour with accuracy proportional
                // to the contour perimeter
                result = cvApproxPoly( contours, sizeof(CvContour), storage,
                    CV_POLY_APPROX_DP, cvContourPerimeter(contours)*0.02, 0 );
                // square contours should have 4 vertices after approximation
                // relatively large area (to filter out noisy contours)
                // and be convex.
                // Note: absolute value of an area is used because
                // area may be positive or negative - in accordance with the
                // contour orientation
                
                if( result->total == 4 &&
                    fabs(cvContourArea(result,CV_WHOLE_SEQ)) > 1000 &&
                    cvCheckContourConvexity(result) )
                {
                    s = 0;
                    
                    for( i = 0; i < 5; i++ )
                    {
                        // find minimum angle between joint
                        // edges (maximum of cosine)
                        if( i >= 2 )
                        {
                            t = fabs(angle(
                            (CvPoint*)cvGetSeqElem( result, i ),
                            (CvPoint*)cvGetSeqElem( result, i-2 ),
                            (CvPoint*)cvGetSeqElem( result, i-1 )));
                            //printf("angulos s %f t %f\n", s,  t );
                            s = s > t ? s : t;
                        }
                    }
                    
                    // if cosines of all angles are small
                    // (all angles are ~90 degree) then write quandrange
                    // vertices to resultant sequence 
                    if( s < 0.3 )
                        for( i = 0; i < 4; i++ )
                            cvSeqPush( squares,(CvPoint*)cvGetSeqElem( result, i ));
                }
                
                // take the next contour
                contours = contours->h_next;
            }
        }
    }
    //cvShowImage( "d",tgray);
    // release all the temporary images
    cvReleaseImage( &gray );
    cvReleaseImage( &pyr );
    cvReleaseImage( &tgray );
    cvReleaseImage( &timg );
    return squares;
}


// the function draws all the squares in the image
void drawSquares( IplImage* img, CvSeq* squares )
{
    CvSeqReader reader;
    IplImage* cpy = cvCloneImage( img );
    int i;
    // initialize reader of the sequence
    cvStartReadSeq( squares, &reader, 0 );
    // read 4 sequence elements at a time (all vertices of a square)
    for( i = 0; i < squares->total; i += 4 )
    {
        CvPoint pt[4], *rect = pt;
        int count = 4;
        // read 4 vertices
        CV_READ_SEQ_ELEM( pt[0], reader );
        CV_READ_SEQ_ELEM( pt[1], reader );
        CV_READ_SEQ_ELEM( pt[2], reader );
        CV_READ_SEQ_ELEM( pt[3], reader );
        //printf("paso x libre %s \n", *rect );
        // draw the square as a closed polyline 
        if (pt[0].x != 1){
         cvPolyLine( cpy, &rect, &count, 1, 1, CV_RGB(0,255,0), 3, CV_AA, 0 );
         //mandar a escribir en GL tambien ahora, k se accede
         amers[0] = pt[0].x; printf(" %d ", (int)amers[0] );
         amers[1] = pt[0].y; printf(" %d ", (int)amers[1] );
         amers[2] = pt[1].x; printf(" %d ", (int)amers[2] );
         amers[3] = pt[1].y; printf(" %d ", (int)amers[3] );
         amers[4] = pt[2].x; printf(" %d ", (int)amers[4] );
         amers[5] = pt[2].y; printf(" %d ", (int)amers[5] );
         amers[6] = pt[3].x; printf(" %d ", (int)amers[6] );
         amers[7] = pt[3].y; printf(" %d \n", (int)amers[7]  );
         //{1,2,3,4,5,6,7};   
         //amers[0] = pt[0]->y;
         calcul_homographie(H_image_vers_modele, H_modele_vers_image, H_GL);
         for(int k = 0 ; k<8;k++) {  printf(" %d ", (int)modele[k] );}printf(" \n ");
         for(int k = 0 ; k<16;k++){  printf(" %d ", (int)H_GL[k] );}printf(" \n ");
         //printf(" %s \n", H_GL  );
        }
    }
    
    // show the resultant image
    cvShowImage( wndname, cpy );
    cvReleaseImage( &cpy );
}

void initRendering() {
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
}


void salirCV (){
     cvReleaseMat(&H_image_vers_modele);
	cvReleaseMat(&H_modele_vers_image);
     cvDestroyWindow( wndname );
    cvDestroyWindow( "d" );
    cvReleaseCapture( &capture );
}

void salir(){
     salirCV();
}

     
static void key(unsigned char key, int x, int y)
{
    switch (key) 
    {
        case 27 : 
        case 'q':
            salir();
            exit(0);
            break;
    }
    glutPostRedisplay();
}

static void display(void){
    GLfloat ambientLight[] = {0.3f, 0.3f, 0.3f, 1.0f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
	
	GLfloat lightColor[] = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat lightPos[] = {-2 * CAJA_LARGO, CAJA_LARGO, 4 * CAJA_LARGO, 1.0f};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0,  ancho, alto);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(45.0, (float)ancho / (float)alto, 1.0, 200.0);
         glOrtho(0,ancho, 0, alto, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	 //glLoadIdentity();
         glLoadMatrixf(H_GL);
	 //glTranslatef(0.0f, 0.0f, -20.0f);
	glColor3f(1., 0., 0.);	
	glBegin(GL_LINE_LOOP);
		glVertex2fv(&modele[0]);
		glVertex2fv(&modele[2]);
		glVertex2fv(&modele[4]);
		glVertex2fv(&modele[6]);
	glEnd();

	glBegin(GL_QUADS);

         //Front face //bleu
	glNormal3f(0.0, 0.0f, 1.0f);
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
    glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
    
    	//Back face
	glNormal3f(0.4, 0.3f, 1.0f);
	glColor3f(0.0f, 0.5f, 0.5f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO );
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO );
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO );
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO );
    	     //Top face //jaune
	glColor3f(1.0f, 1.0f, 0.0f);
	glNormal3f(0.0, 1.0f, 0.0f);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO );
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO );
	glEnd();
/*
	glBegin(GL_QUADS);
	
	//Top face //jaune
	glColor3f(1.0f, 1.0f, 0.0f);
	glNormal3f(0.0, 1.0f, 0.0f);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	
	//Bottom face
	glColor3f(1.0f, 1.0f, 1.0f);
	glNormal3f(0.0, 0.0f, 0.0f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	
	//Left face //vert
	glNormal3f(0.0, 0.0f, 0.0f);
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	
	//Right face //rouge
	glNormal3f(1.0, 0.0f, 0.0f);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	
	//Front face //bleu
	glNormal3f(0.0, 0.0f, 1.0f);
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
    glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, CAJA_ALTO / 2);
    
	//Back face
	glNormal3f(0.4, 0.3f, 1.0f);
	glColor3f(0.0f, 0.5f, 0.5f);
	glVertex3f(-CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(-CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glVertex3f(CAJA_LARGO / 2, -CAJA_ANCHO / 2, -CAJA_ALTO / 2);
	glEnd();
*/
	glutSwapBuffers();
}

static void idle(void) {
       int c;
    //printf("paso x libre s \n" );
       //for( ;;/*i = 0; names[i] != 0; i++*/ )    {
         IplImage* frame=cvQueryFrame(capture);
         if(!frame) return;//break;
        // load i-th image
        //cvCopy( frame, img, 0 );
        
        //img0 = cvLoadImage( names[i], 1 );
        if( !img0 )
        {
            img0 = cvCreateImage( cvGetSize(frame), 8, 3 );
            img0->origin = frame->origin;
            //printf("Couldn't load %s\n", names[i] );
            //break;//continue;
        }
        cvCopy( frame, img0, 0 );
        
        //cvShowImage( "d",img0);
        img = cvCloneImage( img0 );
        cvFlip(img,NULL,1); //OJO SISTEMAS WIN32 VOLTEAN LA IMAGEN /*ATENCION*/ 
        // create window and a trackbar (slider) with parent "image" and set callback
        // (the slider regulates upper threshold, passed to Canny edge detector) 
        cvNamedWindow( wndname, 1 );
        
        // find and draw the squares
        drawSquares( img, findSquares4( img, storage ) );
        
        // wait for key.
        // Also the function cvWaitKey takes care of event processing
        c = cvWaitKey(400);
        // release both images
        cvReleaseImage( &img );
        cvReleaseImage( &img0 );
        // clear memory storage - reset free space position
        cvClearMemStorage( storage );
        if( (char)c == 27 )
            return;//break;
    //}
    glutPostRedisplay();
}

int main(int argc, char** argv)
{
    int i, c;
    // create memory storage that will contain all the dynamic data
    storage = cvCreateMemStorage(0);
    
    capture = cvCaptureFromCAM( argc == 2 ? argv[1][0] - '0' : 0 );
    if( !capture ) return -1;
     ancho =(int)cvGetCaptureProperty( capture,CV_CAP_PROP_FRAME_WIDTH);
     alto  =(int)cvGetCaptureProperty( capture,CV_CAP_PROP_FRAME_HEIGHT);
    cvNamedWindow( "d", 1 );

    glutInit(&argc, argv);
	glutInitWindowPosition(20, 50);
	glutInitWindowSize(ancho, alto);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	
	glutCreateWindow("MODELO");
	initRendering();
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(key);
	
	H_image_vers_modele = cvCreateMat(3, 3, CV_32F);
    
	H_modele_vers_image = cvCreateMat(3, 3, CV_32F);
	
    glutMainLoop();
    salir();
    
    return 0;
}
