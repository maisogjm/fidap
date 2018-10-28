/*      This program displays functional EPI of the form 64x64 matrix
        using X11 windows in multi-slice fashion. Keys and subWindow are added.
        The ar is malloc().
        Image format is 16 bit (short int) array for now.
        Extras: image marker, continuous X,Y position.
        X window skeletons are based on Andre Jesmanowicz's FD.c at MCW.
        Current version uses concatenated multi-slice format.
        Allen W. Song, Lab of Brain and Cognition, NIH, 10.1996 v.1.1.
        Effort is also underway to display volumes with 128 and 256
        matrix in v.1.2.  */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef   SYSV
#include <sys/fcntl.h>
#endif

#include <X11/Xos.h> 
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#include <sys/ioctl.h>

#define IM_HEIGHT 256 
#define IM_ARR    (IM_HEIGHT*IM_HEIGHT)
#define MCOLORS   410                 /* maximum # of colors */
#define NCOLORS   50                 /* default # of colors */
#define N_SPCTR   240                 /* def degree of color spectrum */
#define M_SP_COL  360                 /* max degree of color spectrum */
#define BELT_W    5                  /* reference color belt width */
#define BELT_S    1                   /* color belt sides width */
#define BELT_A    (BELT_W*IM_HEIGHT)
#define NF_MAX    1050                /* Max # of files */
#define STR_L     256                 /* Max length of string */

#define COL_MIN    0
#define MAX_WIDTH  2048  /* Never less then screen max_width */

#define EPX1      64                   /* List of supported EPI image sizes */
#define EPY1      64                   /* Each triple is: */
#define EPS1      (2*EPX1*EPY1)        /* xsize, ysize, filesize */
#define EPX2      128
#define EPY2      128
#define EPS2      (2*EPX2*EPY2)

#define H_SIZE    (IM_ARR)      /* 256x256 image with header */
#define IM_SIZE   (2*H_SIZE)

#define GX_MAX    768                  /* Horizontal size of graph window */
#define GY_MAX    384                  /* Vertical size of graph window */
#define GR_DLX    3                    /* Horizontal delta to right edge */
#define GT_DLY    21                   /* Vertical delta to top edge */
#define GL_DLX    50                   /* Horizontal delta to left edge */
#define GB_DLY    50                   /* Vertical delta to bottom edge */
#define MAT_MAX   15                   /* Maximum array size of graphs */
#define GRID_NUM  5                    /* Maximum grid index */
#define COL_NUM   5                    /* Number of colors */

#define min_max_col(a) ((a) < (256) ? (256) : ((a) > (65280) ? (65280) : (a)))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

/* keys (small subwindows) stuff ----- vvvvvvvvv ------  AJ */

/* 
 #define KFONT     "PellucidaTypewriter12"
 #define TFONT     "PellucidaTypewriter10"
#define KFONT     "*iris-medium-r-normal--12-120-72-72-m-80--ascii"
#define TFONT     "*iris-medium-r-normal--12-120-72-72-m-80--ascii"  */

#define KFONT     "-adobe-times-bold-r-normal--12-*-75-75-p-*-*-*"
#define TFONT     "-adobe-times-bold-r-normal--12-*-75-75-p-*-*-*"

#define PADDINGW  3
#define PADDINGH  8
#define KEY_1_Y   100


struct _key {
           char   *st;
           int    code;
           int    (*fun)();
           Window wid;
           short  x,y,width,height;
  unsigned long   fore,back;
           void   (*func)();
            };

struct _key *key;

#define N_KEYS 9  
#define LAST_K 4
#define kROT   0
#define kHLP   1
#define kAVR   2
#define kDIF   3
#define kNRM   4
#define kAV1   5
#define kAV2   6
#define kIR1   7
#define kIR2   8

int  Ims_rot(), Im_help(), Im_diff(), Ref_im1(), Ref_im2();
int  Im_Aver(), Im_norm(), Av_im1(), Av_im2();

struct _key xtkeys[N_KEYS] = { {"ROTA",  kROT, Ims_rot},
                               {"HELP", kHLP, Im_help},
                               {"AVIM", kAVR, Im_Aver},
                               {"a - b", kDIF, Im_diff},
                               {"Norm", kNRM, Im_norm},
                               {"Average: from this image", kAV1, Av_im1},
                               {"Average: to this one.   ", kAV2, Av_im2},
                               {"Ref. level first image", kIR1, Ref_im1},
                               {"Ref. level last  image", kIR1, Ref_im2} };

char          *kfont, *ffc, *fbc;
unsigned long ForeColor, BackColor, FKeyFore, FKeyBack;
XColor        xcsd, xced;
XFontStruct   *kfontinfo;
Font          keyfont;
int           keywide[N_KEYS], minwide, keyhigh;
GC            Fkeyigc, Fkeygc;
int           invkey    = -1;
int           exp_done[N_KEYS+3];   /* GWindow + subWindow +topWindow + keys */
short   int   a_rot[H_SIZE], rot_nr = 0, rot_direct;
int           sub_W_x, sub_W_y, top_W_x, top_W_y;
int           v_point_x = -20, v_point_y = 0, move_Vpointer = 0;
int           diff_im = 0, Im_1, Im_2, im1_done = 0, extra_im = 0;
int           avr_grp = 0, Av_1, Av_2, av1_done = 0, Av_length = 1;
int           fim_dif = 0, fim_avr = 0, redraw;
int           txtW_ON = 0;


/* keys (small subwindows) stuff ----- ^^^^^^^^^ ------   */

extern double   strtod();
char		*malloc();

XImage          *Load_Any_Arr();
int             STD_colors();

void x_events_loop();
void Resample();
void Put_image();
void plot_line();
void draw_marker();
void scale_up();
void scale_down();
void mat_up();
void mat_down();
void init_mat();
void grid_up();
void grid_down();
void print_plot();
void window_plane();
void graphic_store();
void plotx();
void plx_txt();
void plx_TXT();
void subW_TXT();
void line_color();
void txt_color();
void DrawSubWindow();
void DrawTopWindow();

int    st_4[] = { 1,  2,  3,  4,  3,  2,  1,  0,      /* resampling array */
                  2,  4,  6,  8,  6,  4,  2,  0,
                  3,  6,  9, 12,  9,  6,  3,  0,
                  4,  8, 12, 16, 12,  8,  4,  0,
                  3,  6,  9, 12,  9,  6,  3,  0,
                  2,  4,  6,  8,  6,  4,  2,  0,
                  1,  2,  3,  4,  3,  2,  1,  0,
                  0,  0,  0,  0,  0,  0,  0,  0 };
int    st_2[] = { 1, 2, 1, 0,
                  2, 4, 2, 0,
                  1, 2, 1, 0,
                  0, 0, 0, 0 };

struct S_16_c {
              int red;
              int green;
              int blue;
              };
struct S_16_c Solid_color[16] = {
   {0xffff, 0xffff, 0xffff}, {0xffff, 0x0000, 0x0000},   /* white,  red     */
   {0x0000, 0xffff, 0x0000}, {0x0000, 0x0000, 0xffff},   /* green,  blue    */
   {0x0000, 0xffff, 0xffff}, {0xffff, 0x0000, 0xd3d3},   /* cyan,   magenta */
   {0xffff, 0xffff, 0x0000}, {0xffff, 0x8a8a, 0xffff},   /* yellow, brown   */
   {0x9f9f, 0xd3d3, 0x0000}, {0x0000, 0xffff, 0x9f9f},   /* green,  green   */
   {0x0000, 0x8a8a, 0xffff}, {0x9494, 0x0000, 0xd3d3},   /* blue,   violet  */
   {0xffff, 0x0000, 0x9494}, {0x6969, 0x6969, 0x6969},   /* pink,   gray    */
   {0xaeae, 0xaeae, 0xaeae}, {0x0000, 0x0000, 0x0000}    /* gray,   black   */
  };

short int       color_x11[16] = { -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1 };

XPoint          sm_cir[12] = { {-1,-2},{ 0,-2}, { 1,-2},  /* small circle */
                               { 2,-1},{ 2, 0}, { 2, 1},
                               { 1, 2},{ 0, 2}, {-1, 2},
                               {-2, 1},{-2, 0}, {-2,-1} };
char            *f_name[NF_MAX];
char            *ProgramName = NULL;
char            *display, *geom, **ptr;

int             FIRST = 1, SQUE_NR;
int             x00 = 0, y00 = 0;
int             back_to_main = 0;
XEvent          event, w_event, first_event;

Display         *theDisp;
int             theScreen;
Window          rootW, theWindow, GWindow, subWindow, topWindow;
GC              theGC, txtGC;
unsigned long   fcol, bcol;
Font            mfont;
XFontStruct     *afinfo, *mfinfo;
Visual          *theVisual;
XImage          *theImage, *expImage, *theBelt, *expBelt;
int	        eWIDE, eHIGH, aWIDE, eW, eH, iWIDE, iHIGH;

unsigned short  tmp1[MCOLORS], tmp2[MCOLORS], tmp3[MCOLORS];
int             No_Color_init = 1, No_Grey_init = 1;
int		Dispcells, Planes;
int		repeat_status = 0, ncolors, YES_color = 0, INgrey[MCOLORS];
XColor		MYcol[MCOLORS], MYgrey[MCOLORS], any_col, rgb_col;
Colormap	CMap;
unsigned long	plane_masks[1], pixels[MCOLORS];
unsigned int	nplmsk = 0, ucolors, uucolors;
Pixmap          pxWind;

int             spectrum = N_SPCTR; /* up to 360 degree of colors */
int             Im_Nr = 0, N_im = 0, I_lck = 0, tmp_Nr = 0;
int             N_lck = 0,  old_grid, old_im = -1;
int             x_mag, ref_ar[H_SIZE], av_ar[H_SIZE];
short int       imx[IM_ARR], tmp_imx[IM_ARR], tmp_ar[H_SIZE];
short int       belt_ar[BELT_A], belt_arr[BELT_A];
short int       *ar, *zzz;                                   /* main array */
int             Mltx[MAX_WIDTH], Mlty[MAX_WIDTH];
XPoint          a_line[NF_MAX];

float           del1, coef1, auto_scale;
int             min1;
int             fsize, isize, ar_size;
int             NC;
char            T_name[STR_L];
int             B1_action = 0, Im_frst = 0, fr_index = -2;

int             Argc;
char            **Argv = NULL;

/* these are my own ... EW */

int     plot[MAT_MAX][MAT_MAX][NF_MAX],val[MAT_MAX][MAT_MAX][NF_MAX];
int     pmin[MAT_MAX][MAT_MAX],pmax[MAT_MAX][MAT_MAX],mat,xc,yc,mark;
int     xorigin[MAT_MAX][MAT_MAX],yorigin[MAT_MAX][MAT_MAX];
int     i,j,xpoint,ypoint,npoints,iscale,gx,gy,jimpoint;
int     xspace,yspace,grid_ar[GRID_NUM],grid_index, color_index;
char    plotbuf[10*NF_MAX],pfname[80],fnum[80];
int     xpoint_0 = -1, ypoint_0 = -1; /* old xpoint & ypoint */

int           mytxt, idx, idy;
int           mdx1, mdy1;
char          strp[STR_L];
char          *color[5] = {"yellow","red","cyan","green","cornflowerblue"};

int           im_size, offs;
float         gamm;


/* ---------------- */
   main(argc, argv)
   int  argc;
   char *argv[];
/* ---------------- */
{
   int       i, j, k, m, max1;

   Argc = argc;
   Argv = argv;
   NC = NCOLORS;
   jimpoint = 1;
   ProgramName = argv[0];
   if (argc <  2)  { Syntax(); The_Help(1); }
   display = NULL;
   npoints = 0;
   gamm = 2.;
  
   get_line_args(argc, argv);       /* get command line arguments and files */
  
     im_size = 256;
     x_mag = 1;  /* 256x256 pixels format */

   init_const();

   min1 =  32767;    /* find deviation fot central frame and set scale */
   max1 = -32768;
   i = ypoint*im_size + xpoint;
   for (k=0; k < npoints; k++) {
      m = k * ar_size;
      if ( ar[m+i] < min1 ) min1 = ar[m+i];
      if ( ar[m+i] > max1 ) max1 = ar[m+i];
   }
   del1 = max1 - min1;
   if ( (int) del1 > 1) {
      k = (int) (log(del1*auto_scale)/log(2.) + .5) - iscale; 
      for (i=0; i < abs(k); i++) {       /* scale graph */
         if ( k < 0 ) scale_up();
         else         scale_down();
      }
   }

   min1 =  32767;
   max1 = -32768;
   for (k=0; k < npoints; k++) {    /* find global min and max for main set */
      m = k * ar_size;
      for (i=0; i < ar_size; i++) {
         if ( ar[m+i] < min1 ) min1 = ar[m+i];
         if ( ar[m+i] > max1 ) max1 = ar[m+i];
      }
   }
   del1 = max1 - min1;
   if (del1 < 1.) del1 = 1.;
   coef1 = ( (float) NC - 1.) / del1;


   make_belt(belt_ar, BELT_W, IM_HEIGHT);     /* make reference color belt */
   load_rect_str(belt_ar, belt_arr, BELT_W, IM_HEIGHT,
                 BELT_S, 15, BELT_S, 15, 0, NC);

   main_FD_EPI();          /* Make image in the window and return Drawable */
   
   /* Initialize extra E.C.W. graph stuff here.   */

   window_plane();  /* create graphic window with extra subwindows and keys */
   redraw_graph();                       /* draw frame and text in it */
   DrawSubWindow();
   DrawTopWindow();

   x_events_loop();                              /* Main events check loop */

   exit(0);
 
}       /* --- End of main() --- */

/* -------- */
   Syntax()
/* -------- */
{
  fprintf (stderr, "\n EPI volume display in X11 windows");
  fprintf (stderr, "\n\n Usage: %s [options] vol* or vol1, [vol2, ..., vol%d]\n", ProgramName, NF_MAX);
  fprintf (stderr, "\n Where options are:");
  fprintf (stderr, "\n    -point #_of_point    - initial setting for how large the points in the time course should be");
  fprintf (stderr, "\n    -ngl #_of_gray_levels  - initial number of gray levels [2-%d] (def %d)", MCOLORS, NC);
  fprintf (stderr, "\n    -gam gamma       - gamma correction (1 for no correction)");
  fprintf (stderr, "\n");
}

/* ------------ */
   The_Help(ex)
   int ex;
/* ------------ */
{
  fprintf (stderr, "\n Mouse and Key Events:\n");
  fprintf (stderr,"\n  Program quit      : <q> or <Q>");
  fprintf (stderr,"\n  Change to colors  : <c>");
  fprintf (stderr,"\n  Change to B & W   : <b>");
  fprintf (stderr,"\n  Swap colors       : <s>");
  fprintf (stderr,"\n  Restore original  : right mouse button at image center ");
  fprintf (stderr,"\n  Contrast    : middle or right button - right or left side of the image");
  fprintf (stderr,"\n  Brightness  : middle or right button - top or bottom of the image");
  fprintf (stderr,"\n  First image       : 1");
  fprintf (stderr,"\n  Last image        : l");
  fprintf (stderr,"\n  Next     image    : >");
  fprintf (stderr,"\n  Previous image    : <");
  fprintf (stderr,"\n                      dragging red pointer works too");
  fprintf (stderr,"\n  Scale Plot up         : +");
  fprintf (stderr,"\n  Scale Plot down       : -");
  fprintf (stderr,"\n  Increase Grid Spacing : G");
  fprintf (stderr,"\n  Decrease Grid Spacing : g");
  fprintf (stderr,"\n  Matrix size up  : u");
  fprintf (stderr,"\n  Matrix size down : d");
  fprintf (stderr,"\n  Exact matrix size     : N #of_size <CR> (1 to %d only)", MAT_MAX);
  fprintf (stderr,"\n  Save minigraph in ASCII file   : press <p>");
  fprintf (stderr,"\n  Save current image to a file   : press <S>");
  fprintf (stderr,"\n  Position frame in the image    : press left button in the image area.");
  fprintf (stderr,"\n  Center frame on desired pixel  : press left button on desired minigraph.");
  fprintf (stderr,"\n  Change to differential display : press [a - b] key. Set first and last image.");
  fprintf (stderr,"\n\n ");
  if (ex == 1) exit(1);
}

/* ------------------------- */
   get_line_args(argc, argv)
   int  argc;
   char *argv[];
/* ------------------------- */
{
   register int i, j, k, nopt, nnn;
   int          sp;
   float        fff;

   nopt = 0;
   for (i = 1; i < argc; i++) { /* ------- Options ------- */

      if (!strncmp(argv[i], "-d", 2)) {          /* display */
         if (++i >= argc) { Syntax (); exit(2); }
         display = argv[i];
         nopt++; nopt++;
         continue;
      }
      if (!strncmp(argv[i], "-point", 3)) {          /* initial geometry */
         if (++i >= argc) { Syntax (); exit(2); }
         jimpoint = atoi(argv[i]);
         nopt++; nopt++;
         continue;
      }
      if (!strncmp(argv[i], "-h", 2)) {          /* help */
         Syntax();
         The_Help(1);
      }
      if (!strncmp(argv[i], "-ngl", 3)) {         /* # of colors */
         if (++i >= argc) { Syntax (); exit(2); }
         nnn = atoi(argv [i]);
         if (nnn < 2 || nnn > MCOLORS ) { Syntax(); exit(2); }
         NC = nnn;
         nopt++; nopt++;
         continue;
      }
      if (strncmp(argv [i], "-sp", 3) == 0) {
         if (++i >= argc) { Syntax(); exit(2); }
         ptr = argv;
         sp = strtod(argv[i], ptr);
         if ( **ptr || (sp < 1) || (sp > M_SP_COL) ) {
            fprintf (stderr, "\n !!! Wrong degree value: -sp %d !!!\n",sp);
            fprintf (stderr, "\n !!! For help type: %s -h \n\n", ProgramName);
            exit(3);
         }
         spectrum = sp;
         nopt++; nopt++;
         continue;
      }
      if (strncmp(argv [i], "-gam", 4) == 0) {
         if (++i >= argc) Syntax ();
         ptr = argv;
         fff = strtod(argv[i], ptr);
         if ( **ptr || (fff <= 0.) ) {
            fprintf (stderr, "\n !!! Wrong gamma value: %g !!!\n\n",fff);
            Syntax();
         }
         gamm = fff;
         nopt++; nopt++;
         continue;
      }
      if (strncmp(argv [i], "-num", 4) == 0) {
         if (++i >= argc) { Syntax(); exit(2); }
         ptr = argv;
         npoints = strtod(argv[i], ptr) + .5;
         if ( **ptr || (npoints < 1)) {
            fprintf (stderr, "\n !!! Too few images specified !!!\n\n");
            exit(1); /* now symbolic for min npoints = 1 .  */
         }
         nopt++; nopt++;
         continue;
      }
      if (strncmp(argv [i], "-im1", 4) == 0) {
         if (++i >= argc) { Syntax(); exit(2); }
         ptr = argv;
         Im_frst = strtod(argv[i], ptr) + .5;
         if ( **ptr || (Im_frst < 1)) {
            fprintf (stderr, "\n !!! First_image_# < 1 in -im1 !!!\n\n");
            exit(1);
         }
         nopt++; nopt++;
         continue;
      }
      if (strncmp(argv [i], "-mag", 4) == 0) {
         fprintf (stderr, "\n !!! No more -mag option in use !!!");
         fprintf (stderr, "\n !!! Resize window using mouse  !!!\n\n");
         if (++i >= argc) { Syntax(); exit(2); }
         nopt++; nopt++;
         continue;
      }

   }
   nopt++;                                   /* Files to read (minimum one) */

   if ( nopt > (argc-1) || nopt < (argc - NF_MAX) ) {  /* Nr of files check */
      fprintf (stderr, "\n Wrong # of files. %d files entered :\n", argc-nopt);
      for(i=nopt, j=1; i < argc; i++, j++)
         fprintf (stderr, "  %3d -  %s\n", j, argv[i]);
         Syntax(); exit(2);
   }

   N_im = argc-nopt;                                          /* # of images */
   if ( Im_frst > N_im ) {
      fprintf (stderr, "\n !!! First_im_# in -im1 is bigger then number of images !!!\n\n");
      Syntax(); exit(2);
   }
   for (i=0; i < N_im; i++) f_name[i] = Argv[i+nopt];
   for (i=0; i < Im_frst-1; i++) f_name[i] = f_name[Im_frst-1];

               /* read and check the length of the first file for validity */
   isize = 0;
   printf("\n\n Reading volume %s\n", f_name[0]);
   if (k = read_iqm(f_name[0], &isize, tmp_ar)) {
      fprintf (stderr, "\n Problem with file: %s\n", f_name[0]);
      Syntax(); exit(2);
   }
 /*  if ( (isize != EPS1) && (isize != EPS2) && (isize != IM_SIZE) ) {
      fprintf (stderr, "\n\n !!! File %s has wrong format !!!\n", f_name[0]);
      Syntax(); exit(2);
    }   */         /* disabled for arbitrary volume display   */

   fsize = isize;
   if (fsize == IM_SIZE) {
      ar_size = IM_ARR;    /* 256x256 + header image */
      offs = 0;
   }
   else {
      ar_size = fsize/2;   /* AxA and no header */
      offs    = 0;
   }

   ar = (short int *) malloc((unsigned) ((ar_size*N_im)*sizeof(short int)));
   zzz = (short int *) malloc((unsigned) ((IM_ARR*N_im)*sizeof(short int)));

   for (i=0; i < ar_size; i++) ar[i] = tmp_ar[offs+i]; /* reload first image */

   for (i=1; i < N_im; i++) {                            /* Read image files */
      printf(" Reading volume %s\n", f_name[i]);
      isize = 0;
      if ( k = read_iqm(f_name[i], &isize, tmp_ar) ) {
         fprintf (stderr, "\n Problem [%d] with file: %s\n", k, f_name[i]);
         Syntax(); exit(2);
      }
      if ( isize != fsize) {  /* check lengthes of other files */ 
         fprintf (stderr, "\n\n !!! File %s has different format !!!\n",
                  f_name[i]);
         Syntax(); exit(2);
      }  
      k = i*ar_size;
      for (j=0; j < ar_size; j++) ar[k+j] = tmp_ar[offs+j]; /* reload data */
   }
   if (!npoints || npoints > N_im) npoints = N_im;
}

/* ------------------------------------------- Main link from main() C */
   main_FD_EPI()
/* ------------------------------------------------------------------- */
{
   int        i, z;
   XGCValues  gcv;

   ucolors = SQUE_NR = ncolors = NC;
   uucolors = NC - 1;
   expImage = NULL;

   for (i=0; i<IM_HEIGHT; i++) {   /* init MResize() lookup tables */
      Mltx[i] = i;
      Mlty[i] = i;
   }

   top_nr_name(f_name[0]);

/* ---------- Open the display. --------- */

   if ( (theDisp=XOpenDisplay(display)) == NULL) {
     fprintf(stderr, "%s: Can't open display (are X11 windows running ?). \n"
			,ProgramName);
     exit(1);
   }

   theScreen = DefaultScreen(theDisp);		/* Set Graphic variables */
   rootW     = RootWindow(theDisp,theScreen);
   theGC     = DefaultGC(theDisp,theScreen);
   theVisual = DefaultVisual(theDisp,theScreen);
   CMap      = DefaultColormap(theDisp, theScreen);
   Planes    = DisplayPlanes(theDisp, theScreen);

   if ((Dispcells = DisplayCells(theDisp, theScreen)) <= 200)
      FatalError("This program requires min 8 bits per pixel. ");

   if (!(XAllocNamedColor(theDisp, CMap, "black", &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. ");
   fcol = any_col.pixel;

   if (!(XAllocNamedColor(theDisp, CMap, "grey35", &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. ");
   bcol = any_col.pixel;

   if (!(XAllocColorCells(theDisp, CMap, True, plane_masks, nplmsk,
      pixels, ucolors))) FatalError ("Reduce number of gray levels. AWS ");

   colmap_init();                              /* Initialize color map */

         /* Create theImage from ar[] #1. Do NOT use here Put_image() !  */
   for (z=0; z< N_im; z++) {
   for (i=0; i < ar_size; i++)                   /* renormalize data */
      tmp_ar[i+offs] = ar[z*ar_size + i];
   Resample(&tmp_ar[offs]);
      for (i=0; i< IM_ARR; i++)
      zzz[z*IM_ARR + i] = imx[i];
   }

   for (i=0; i < ar_size; i++)                   /* renormalize data */
      tmp_ar[i+offs] = (int)((float)(ar[i]-min1)*coef1 + .5);
   Resample(&tmp_ar[offs]);

   theImage = Load_Any_Arr(imx,IM_HEIGHT,IM_HEIGHT);

   theBelt  = Load_Any_Arr(belt_arr,BELT_W,IM_HEIGHT);

   eWIDE = theImage->width;	eHIGH = theImage->height;
   eW = theBelt->width;	eH = theBelt->height;
   aWIDE = eWIDE + eW;

   if (!(afinfo = XLoadQueryFont(theDisp, KFONT)))
      if (!(afinfo = XLoadQueryFont(theDisp, "fixed"))) {
         sprintf(strp, "Can't open %s or fixed fonts\n", KFONT);
         FatalError (strp);
      }
   
   mfont=afinfo->fid;
   XSetFont(theDisp,theGC,mfont);
   XSetForeground(theDisp,theGC,fcol);
   XSetBackground(theDisp,theGC,bcol);

   CreateMainWindow(geom, Argc, Argv);	/* Set theWindow properties */

   gcv.function = GXcopy;		/* Extra, Small Text Graphic Context */
   gcv.foreground = fcol;
   txtGC = XCreateGC(theDisp, theWindow, GCForeground|GCFunction, &gcv);
   if (!(mfinfo = XLoadQueryFont(theDisp, TFONT))) /*Set font*/
      if (!(mfinfo = XLoadQueryFont(theDisp, "fixed"))) {
         sprintf(strp, "Can't open %s or fixed fonts\n", TFONT);
         FatalError (strp);
      }
   mfont=mfinfo->fid;
   XSetFont(theDisp,txtGC,mfont);
   XSetForeground(theDisp,txtGC,fcol);
   XSetBackground(theDisp,txtGC,bcol);

   New_Cursor(theDisp, theWindow, XC_left_ptr, "red", "white"); /* cursor */

				/* Adjust size before Show */
   MResize(&expImage,&eWIDE,&eHIGH,theImage,eWIDE,eHIGH);
   MResize(&expBelt,&eW,&eH,theBelt,eW,eH);

   XSelectInput(theDisp, theWindow, ExposureMask | KeyPressMask
           | ButtonPressMask | ButtonReleaseMask | StructureNotifyMask);

   XMapWindow(theDisp,theWindow);		/* Show theImage first time */

   while (1) {
      XNextEvent(theDisp, &first_event);        /* Wait for window on */
      switch (first_event.type) {
         case Expose: {
	   XExposeEvent *e = (XExposeEvent *) &first_event;
           if (e->window == theWindow) {
	     XPutImage(theDisp, theWindow, theGC, expImage,
		 e->x, e->y, e->x, e->y, e->width, e->height);
	     XPutImage(theDisp, theWindow, theGC, expBelt,
	         0, 0, eHIGH, 0, eW, eH);
             Allow_smaller_im(Argc, Argv);
	   }
         return;
	 }

         case KeyPress:
         case ButtonPress:
         case ButtonRelease:
         case ConfigureNotify:
         case CirculateNotify:
         case MapNotify:
         case DestroyNotify:
         case GravityNotify:
         case ReparentNotify:
         case UnmapNotify:
         break;

         default:                /* ignore unexpected events */
         break;
      }
   }
}        /* end of main_FD_EPI() */

/* ----------------------- */
   void x_events_loop()
/* ----------------------- */
{
   register int i, n, m;

   for (i=LAST_K; i<N_KEYS; i++) XUnmapWindow(theDisp,  key[i].wid);
   back_to_main = 0;
   while (1) {
     if (repeat_status) {       /* For repeated action when button pressed */
       n = max(0, (120 - ncolors)*800);
       for(i=0; i < n; i++) m = 5*ncolors;
       if (XCheckWindowEvent(theDisp, theWindow, ButtonReleaseMask, &w_event))
         event = w_event;
     }
     else {
       XNextEvent(theDisp, &event);  /* single event loop */
     }
     HandleEvent(&event);
     if(back_to_main) return;
   }
}

/* ------------------ */
   HandleEvent(event)
   XEvent *event;
/* ------------------ */
{
   Window wind;

 switch (event->type) {

   case Expose: {
      XExposeEvent *e = (XExposeEvent *) event;
      wind = e->window;
      if (wind == theWindow) {
         XPutImage(theDisp, theWindow, theGC, expImage,
                   e->x, e->y, e->x, e->y, e->width, e->height);
         XPutImage(theDisp, theWindow, theGC, expBelt,
                   0, 0, eHIGH, 0, eW, eH);
      }
      else if (wind == subWindow) {
         DrawSubWindow();
      }
      else if (wind == topWindow) {
         DrawTopWindow();
      }
      else {
         int  i;
         for (i=0; i < N_KEYS; i++)
            if (key[i].wid==wind && i > 2) DrawKey(i);
         }
      }
       break;

   case	KeyPress: {
      XKeyEvent *key_event = (XKeyEvent *) event;
      char      buf[128];
      int       i;
      KeySym    ks;
      XComposeStatus status;

      buf[0] = 0;
      XLookupString(key_event, buf, 128, &ks, &status);
      if (buf[0]=='q' || buf[0]=='Q') exit(0);
      if ( (I_lck == 0) && (buf[0]=='I') ) {    /* init image number */
         I_lck = 1;
         New_Cursor(theDisp, theWindow, XC_hand2, "red", "white");
         New_Cursor(theDisp,   GWindow, XC_hand2, "yellow", "red");
         break;
      }
      if ( (N_lck == 0) && (buf[0]=='N') ) {  /* init grid resol. number */
         N_lck = 1;
         old_grid = mat;
         New_Cursor(theDisp, theWindow, XC_hand2, "red", "white");
         New_Cursor(theDisp,   GWindow, XC_hand2, "yellow", "red");
         break;
      }

      if ( (I_lck == 1) && (buf[0] == 13) ) { /* make image number */
         if ( tmp_Nr >= 0 ) {
            int aaa = avr_grp, ddd = diff_im;

            if ( tmp_Nr > 0 ) Im_Nr = min(N_im - 1, tmp_Nr - 1);
            New_Cursor(theDisp, theWindow, XC_left_ptr, "red", "white");
            New_Cursor(theDisp,   GWindow, XC_left_ptr, "blue", "yellow");
            XClearWindow(theDisp, GWindow);
            diff_im = fim_dif;
            avr_grp = fim_avr;
            if ( Im_Nr >= npoints ) {
               avr_grp = 0;
               diff_im = 0;
            }
            else 
               if ( Im_Nr > npoints - Av_length )
                  Im_Nr = npoints - Av_length;
            redraw =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
            Put_image(Im_Nr);
            if ( redraw ) redraw_graph();
            DrawSubWindow();
            DrawTopWindow();
            I_lck = 0;
            tmp_Nr = 0;
         }
         break;
      }
      if ( (N_lck == 1) && (buf[0] == 13) ) { /* make grid resol. number */
         if ( tmp_Nr >= 0 ) {
            if ( tmp_Nr > 0 ) mat = min(MAT_MAX, tmp_Nr);
            New_Cursor(theDisp, theWindow, XC_left_ptr, "red", "white");
            New_Cursor(theDisp,   GWindow, XC_left_ptr, "blue", "yellow");
            if (mat != old_grid) {
               init_mat();
               redraw_graph();
            }
            Put_image(Im_Nr);
            DrawSubWindow();
            N_lck = 0;
            tmp_Nr = 0;
         }
         break;
      }
      if ( (I_lck == 1) || (N_lck == 1) ) {   /* make a number from digits */
         if ( isdigit(buf[0]) ) {
            i = buf[0] - 48;
            tmp_Nr = min(10000, 10*tmp_Nr + i);
         }
         break;
      }
      else {
         switch (buf[0]) {
            case '^': {
               New_Cursor(theDisp, theWindow, XC_left_ptr, "red", "white");
               break;
            }
            case 's': {
               c_swap();
               break;
            }
            case 'S': {
               save_act_im();
               break;
            }
            case 'c': {
               YES_color = 1;
               colmap_init();
               break;
            }
            case 'b': {
               YES_color = 0;
               colmap_init();
               break;
            }
            case '1': {
               int aaa = avr_grp, ddd = diff_im;

               Im_Nr = 0;
               diff_im = fim_dif;
               avr_grp = fim_avr;
               redraw =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
               XClearWindow(theDisp, GWindow);
               Put_image(Im_Nr);
               if ( redraw ) redraw_graph();
               DrawSubWindow();
               DrawTopWindow();
               discard(KeyPressMask, event);
               break;
            }
            case 'l': {
               int max_im = N_im - 1;
               int aaa = avr_grp, ddd = diff_im;

               if ( max_im >= npoints ) {
                  avr_grp = 0;
                  diff_im = 0;
               }
               else {
                  diff_im = fim_dif;
                  avr_grp = fim_avr;
                  max_im = npoints - Av_length;
               }
               redraw =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
               Im_Nr = max_im;
               XClearWindow(theDisp, GWindow);
               Put_image(Im_Nr);
               if ( redraw ) redraw_graph();
               DrawSubWindow();
               DrawTopWindow();
               discard(KeyPressMask, event);
               break;
            }
            case '>': {
               int max_im, i_tmp = Im_Nr;
               int aaa = avr_grp, ddd = diff_im;

               max_im = npoints - Av_length;

               if ( Im_Nr >= npoints ) {
                  if ( Im_Nr < (N_im - 1) ) Im_Nr += 1;
               }
               else if ( Im_Nr < max_im ) Im_Nr += 1;
               else if ( npoints < N_im ) {
                  avr_grp = 0;
                  diff_im = 0;
                  Im_Nr = npoints;
               }
               if ( i_tmp != Im_Nr ) {
                  redraw =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
                  XClearWindow(theDisp, GWindow);
                  Put_image(Im_Nr);
                  if ( redraw ) redraw_graph();
                  DrawSubWindow();
                  DrawTopWindow();
                  discard(KeyPressMask, event);
               }
               break;
            }
            case '<': {
               int min_im = 0, i_tmp = Im_Nr;
               int aaa = avr_grp, ddd = diff_im;

               if ( Im_Nr > npoints ) Im_Nr -= 1;
               else {
                  diff_im = fim_dif;
                  avr_grp = fim_avr;
                  if(Im_Nr > min_im) Im_Nr -= 1;
                  if ( avr_grp )
                     if ( Im_Nr > (npoints - Av_length) )
                        Im_Nr = npoints - Av_length;
               }
               if ( i_tmp != Im_Nr ) {
                  redraw =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
                  XClearWindow(theDisp, GWindow);
                  Put_image(Im_Nr);
                  if ( redraw ) redraw_graph();
                  DrawSubWindow();
                  DrawTopWindow();
                  discard(KeyPressMask, event);
               }
               break;
            }
                   /* + and - keys used for adjusting scale */
            case '-': {
               scale_down();
               redraw_graph();
               DrawSubWindow();
               break;
            }
            case '+': {
               scale_up();
               redraw_graph();
               DrawSubWindow();
               break;
            }
                    /* p used to write plot to file */
            case 'p': {
               print_plot();
               discard(KeyPressMask, event);
               break;
            }
                   /* d and u keys used for adjusting matrix */
            case 'd': {
               mat_down();
               discard(KeyPressMask, event);
               break;
            }
            case 'u': {
               mat_up();
               discard(KeyPressMask, event);
               break;
            }
                   /* g and G keys used for adjusting grid lines */
            case 'g': {
               grid_down();
               discard(KeyPressMask, event);
               break;
            }
            case 'G': {
               grid_up();
               discard(KeyPressMask, event);
               break;
            }
            case 'r': {
               color_index++;
               if (color_index >= COL_NUM) color_index = -1;
               if (color_index < 1) mark = 1 - mark;
               redraw_graph();
               DrawSubWindow();
               discard(KeyPressMask, event);
               break;
            }
            case 'R': {
               fr_index--;
               if (fr_index < -16) fr_index = -1;
               Put_image(Im_Nr);
               discard(KeyPressMask, event);
               break;
            }
         }
      } /* end I_lck = 0 */

   }  /* end of KeyPress */
	break;

   case	ButtonPress: {
      XButtonPressedEvent *b = (XButtonPressedEvent *) event;
      int delta = 1;

      wind = b->window;
      if (wind == theWindow ) {

         if (b->button == Button1 && repeat_status == 0) {
            if (b->x > eWIDE) break;
            B1_action = 1;
            Track_Cursor(b->x, b->y);
            discard(ButtonPressMask, event);
            break;
         }

         if ( b->button == Button2 ) {        /* contrast and brightness */
            repeat_status = 1;
            if (b->x <   eWIDE/3) {
               c_squeeze( delta); break;
            }
            if ((b->x > 2*eWIDE/3) && (b->x < eWIDE)) {
               c_squeeze(-delta); break;
            }
            if ((b->y <   eHIGH/3) && (b->x < eWIDE)) {
               c_bright(-delta); break;
            }
            if ((b->y > 2*eHIGH/3) && (b->x < eWIDE)) {
               c_bright( delta); break;
            }
            if ((b->y >   eHIGH/2) && (b->x > eWIDE) && (b->x < eWIDE + eW)) {
               c_rotate(-delta); break;
            }
            if ((b->y <   eHIGH/2) && (b->x > eWIDE) && (b->x < eWIDE + eW)) {
               c_rotate(delta); break;
            }
         }
   
         if (b->button == Button3 && repeat_status == 0) {
            delta = 1;
            if (b->x <   eWIDE/3) {
               c_squeeze( delta); break;
            }
            if ((b->x > 2*eWIDE/3) && (b->x < eWIDE)) {
               c_squeeze(-delta); break;
            }
            if ((b->y <   eHIGH/3) && (b->x < eWIDE)) {
               c_bright(-delta); break;
            }
            if ((b->y > 2*eHIGH/3) && (b->x < eWIDE)) {
               c_bright( delta); break;
            }
            if ((b->y >   eHIGH/3) && (b->y < 2*eHIGH/3) &&
                (b->x >   eWIDE/3) && (b->x < 2*eWIDE/3)) {
               if (YES_color) No_Color_init = 1;
               else           No_Grey_init  = 1;
               colmap_init(); break;
            }
            if ((b->y >   eHIGH/2) && (b->x > eWIDE) && (b->x < eWIDE + eW)) {
               c_rotate(-delta); break;
            }
            if ((b->y <   eHIGH/2) && (b->x > eWIDE) && (b->x < eWIDE + eW)) {
               c_rotate(delta); break;
            }
         }
      }   /* for theWindow */
      else if (wind == GWindow) {
         int i, j, k1, k2;
      
         if (b->button == Button1) {
            k1 = GX_MAX/mat;
            k2 = GY_MAX/mat;
            i = (b->x - GL_DLX + k1)*mat/GX_MAX;
            j = (b->y - GT_DLY + k2)*mat/GY_MAX;
            xpoint += i - (mat + 2)/2;
            ypoint += j - (mat + 2)/2;
            if (xpoint < 0)        xpoint += im_size;
            if (xpoint >= im_size) xpoint -= im_size;
            if (ypoint < 0)        ypoint += im_size;
            if (ypoint >= im_size) ypoint -= im_size;

            if ( (xpoint == xpoint_0 ) && (ypoint == ypoint_0) ) ;
            else {
               Put_image(Im_Nr);
               redraw_graph();
               DrawSubWindow();
               xpoint_0 = xpoint;
               ypoint_0 = ypoint;
            }
            discard(ButtonPressMask, event);
         }
         break;
      }   /* for GWindow */
      else if (wind == subWindow) {
         if (b->button == Button1) {
            if ( (abs(b->x - v_point_x) < 5 ) &&
                 (abs(b->y - v_point_y) < 9 ) ) {
               move_Vpointer = 1;
               Track_Vpointer();
            }
         }
      }   /* for subWindow */
   
      if ( (b->button == Button1) || (b->button == Button3) ) {
         int  i;
         if      (b->button == Button1) rot_direct = -1;
         else if (b->button == Button3) rot_direct =  1;
         for (i=0; i < N_KEYS; i++) {
            if (key[i].wid==wind && i > 2) { InvertKey(i); invkey=i; }
         }
      }   /* for keys */
   }   /* end of ButtonPress */
   break;

   case	ButtonRelease: {
      XButtonReleasedEvent *b = (XButtonReleasedEvent *) event;

      wind = b->window;

      if (wind == theWindow ) {
         if (B1_action == 1) {
            B1_action = 0;
            xpoint = Mltx[b->x]/x_mag;
            ypoint = Mlty[b->y]/x_mag;

            if ( (xpoint == xpoint_0 ) && (ypoint == ypoint_0) ) ;
            else {
               Put_image(Im_Nr);
               redraw_graph();
               DrawSubWindow();
               xpoint_0 = xpoint;
               ypoint_0 = ypoint;
            }
         }
         if (b->button == Button2) {
            repeat_status = 0;
         }
      }
      else if (wind == subWindow ) {
         discard_Key(wind, ButtonPressMask, event);
         discard_Key(wind, ButtonReleaseMask, event);
      }
      else if ( (b->button == Button1) || (b->button == Button3) ) {
         int i;
         for (i=0; i < N_KEYS; i++) {
            if (wind == key[i].wid && i > 2) LetGoKey(invkey);
            discard_Key(wind, ButtonPressMask, event);
            discard_Key(wind, ButtonReleaseMask, event);
         }
         invkey = -1;
      }
   }
   break;

   case	ConfigureNotify: {
      XConfigureEvent *conf_event = (XConfigureEvent *) event;
      int wb;
      wb = conf_event->height*BELT_W/IM_HEIGHT;
      if ( (conf_event->window == theWindow) && 
           (conf_event->height != eHIGH)) {
         MResize(&expImage, &eWIDE, &eHIGH, theImage,
                  conf_event->height, conf_event->height);
         MResize(&expBelt, &eW, &eH, theBelt,
                  wb, conf_event->height);
      }
   }
   break;


   case CirculateNotify:
   case MapNotify:
   case DestroyNotify:
   case GravityNotify:
   case ReparentNotify:
   case UnmapNotify:
      break;

   default:		/* ignore unexpected events */
      break;
 }   /* end of switch event->type */
}   /* end of HandleEvent() */

/* It loads color strip made by make_belt() & adds unit color stips */
/* ------------------------------------------------------------- */
   load_rect_str(idata, ipx, x, y, i1, x1, i2, x2, amin, nshade)
   short int idata[], ipx[];
   int       x, y, i1, i2, x1, x2, amin, nshade;
/* ------------------------------------------------------------- */
{
   register int i, ii, j, k;
   int          min0, max1;
   float        adel0, coef;

   min0 =  32767;
   max1 = -32768;
   ii = x*y;
   for (i=0; i < ii; i++) {
      if ( idata[i] < min0 ) min0 = idata[i];
      if ( idata[i] > max1 ) max1 = idata[i];
   }
   adel0 = max1 - min0;

   /* set image shades */

   if (adel0 < 1.) adel0 = 1.; 
   coef = ( (float) nshade - 1. ) / adel0;

   k = x - i2;
   for (i=0; i < y; i++) {
      ii = i*x;
      for (j=0; j < i1; j++)
         ipx[ii+j] = -x1;
      for (j=x-i2; j < x; j++)
         ipx[ii+j] = -x2;
      for (j=i1; j < k; j++)
         ipx[ii+j] = (int) ( (float) (idata[ii+j] - min0) * coef + .5) + amin;
   }
}

/* ------------------- */       /* It makes data belt (ramp of indices) */
   make_belt(id, x, y)
   short int *id;
   int       x, y;
/* ------------------- */
{
   register int i, j, ii;

   for(i=0; i < y; i++) {
      ii = i*x;
      for(j=0; j < x; j++)
         id[ii+j] = x - i;
   }
}

/* ------------------------------ */ /* It resizes images universally  */
   MResize(emage,eW,eH,image,w,h)    /* for 12 and 8 bit planes. */
   int     w, h, *eW, *eH;
   XImage  *image,*(*emage);
/* ------------------------------ */
{
   int          lt[MAX_WIDTH];
   register int iy, ex, ey, iW, iH, iG, iN, iE, w2, tmp;
   char         *ximag;
   char         *Ep, *El, *Ip, *Il, *Id;

    iN = image->bits_per_pixel/8;    /* 1 or 2, pointer increment */
    iE = iN - 1;                     /* 0 or 1 for next pointer */
    iG = image->bytes_per_line;      /* Number of bytes in line */
    w2 = w * iN;
    iW = image->width;
    iH = image->height;
    if (w==iW && h==iH) {       /* very special case */
        if (*emage != image) {
            if (*emage) XDestroyImage(*emage);
            *emage = image;
            *eW = iW;  *eH = iH;
        }
    }

    else {                              /* have to do some work */

        /* first, kill the old emage, if one exists */
        if (*emage && *emage != image) {
           free((*emage)->data);  (*emage)->data = NULL;
           XDestroyImage(*emage);
        }

        /* create emage of the appropriate size */

        *eW = w;  *eH = h;
        ximag = (char *) malloc(w2 * h);

        *emage = XCreateImage(theDisp, theVisual, Planes, ZPixmap, 0, ximag,
                        w, h, 8, w2);

        if (!ximag || !*emage) {
           fprintf(stderr,"ERROR: unable to create a %dx%d image\n",w,h);
           exit(1);
        }

        El = Ep = (char *) (*emage)->data;
        Id = (char *) image->data;

        /* ---- Set look up table and Reframe fast ---- */

        for (ex = 0; ex < w; ex++)  {
            tmp = (iW * ex)/w;
            lt[ex] = iN * tmp;
            Mltx[ex] = tmp;
        }

        for (ey = 0; ey < h; ey++, El += w2) {
            iy = (iH * ey)/h;
            Mlty[ey] = iy;
            Ep = El;
            Il = Id + iG * iy;
            for(ex = 0; ex < w; ex++, Ep += iN) {
                Ip = Il + lt[ex];
                *Ep = *Ip;
                *(Ep+iE) = *(Ip+iE);
            }
        }
    }
}

/* ---------------------------- */
   Allow_smaller_im(argc, argv)
   int  argc;
   char *argv[];
/* ---------------------------- */
{
   XSizeHints           hints;

   hints.min_height = 0;
   hints.min_width = 0;
   hints.flags = PMinSize;

   hints.min_aspect.x = IM_HEIGHT + BELT_W;
   hints.max_aspect.x = IM_HEIGHT + BELT_W;
   hints.min_aspect.y = IM_HEIGHT;
   hints.max_aspect.y = IM_HEIGHT;
   hints.flags |= PAspect;

   hints.max_height = DisplayHeight(theDisp, theScreen);
   hints.max_width = hints.max_height *
		 (IM_HEIGHT + BELT_W)/IM_HEIGHT;
   hints.flags |= PMaxSize;
   XSetStandardProperties(theDisp, theWindow, T_name, T_name, None,
			  argv, argc, &hints);
}

/* ---------------------------------- */
   CreateMainWindow(geom, argc, argv)
   char *geom;
   int  argc;
   char *argv[];
/* ---------------------------------- */
{
   XClassHint		class;
   XSetWindowAttributes attr;
   unsigned int		attrmask;
   XSizeHints		hints;
   int 			mask, x, y, w, h, i;

   x = y = w = h = 1;
   mask = XParseGeometry(geom, &x, &y, &w, &h);

   if (mask & HeightValue) { 
     i = IM_HEIGHT + BELT_W;
     if ( i > MAX_WIDTH ) h = (MAX_WIDTH * IM_HEIGHT) / i;
     eH = eHIGH = h;
     eWIDE = eHIGH;
     eW = eH * BELT_W / IM_HEIGHT;
     aWIDE = eWIDE + eWIDE*BELT_W/IM_HEIGHT;
   }

   if (mask & XValue || mask & YValue) hints.flags = USPosition;  
   else hints.flags = PPosition;

   hints.flags |= USSize;

   if (mask & XValue && mask & XNegative) 
        x = XDisplayWidth(theDisp, theScreen)-eWIDE-abs(x);
   if (mask & YValue && mask & YNegative) 
        y = XDisplayHeight(theDisp, theScreen)-eHIGH-abs(y);

   class.res_name  = "ImageX";
   class.res_class = "ImageX";

   hints.min_aspect.x = IM_HEIGHT + BELT_W;
   hints.max_aspect.x = IM_HEIGHT + BELT_W;
   hints.min_aspect.y = IM_HEIGHT;
   hints.max_aspect.y = IM_HEIGHT;
   hints.flags |= PAspect;

   hints.x = x;			hints.y = y;
   hints.width = aWIDE;		hints.height = eHIGH;
   hints.max_width = aWIDE;	hints.max_height = eHIGH;
   hints.flags |= PMaxSize;

   hints.min_width = aWIDE;	hints.min_height = eHIGH;;
   hints.flags |= PMinSize;

   attr.background_pixel = bcol;
   attr.border_pixel     = fcol;
   attrmask = CWBackPixel | CWBorderPixel;

   theWindow = XCreateWindow(theDisp, rootW, x, y, aWIDE, eHIGH, 2,
	CopyFromParent, CopyFromParent, CopyFromParent, attrmask, &attr);

   if (!theWindow)
      FatalError("unable to open window (X11 windows not running).");

   XSetClassHint(theDisp, theWindow, &class);
   XSetStandardProperties(theDisp, theWindow, T_name, T_name, None,
			  argv, argc, &hints);
}

/* ----------------- */	/* Set color or grey array of the Image	*/
   colmap_init()
/* ----------------- */
{
   if (YES_color)  color_init(COL_MIN, spectrum);
   else            grey_init();
}

/* -------------------- */       /* Set grey array of the Image  */
   grey_init()
/* -------------------- */
{
   register int i, k, m;
   float a;

   a = 200./ncolors;

   if (No_Grey_init) {
      for (i=0; i < ncolors; i++) {
         k = 255.*pow((a*i+55.)/255., gamm) + .5;
         m = min_max_col(k << 8);
         INgrey[i] = m;
         MYgrey[i].pixel = pixels[i];
         MYgrey[i].red   = m;
         MYgrey[i].green = m;
         MYgrey[i].blue  = m;
         MYgrey[i].flags = DoRed|DoGreen|DoBlue;
      }
      No_Grey_init = 0;
   }
   XStoreColors(theDisp, CMap, MYgrey, ncolors);
}

/* ----------------- */	/* Set color array of the Image	*/
   color_init(a1,a2)    /* OK but sharp edge on green near yellow */
   int	a1, a2;
/* ----------------- */
{
   double da, an, c, n, s, sb, cb, ak, ab;
   register int k, i, m, r, g, b;

   if (No_Color_init) {
      FIRST = 1;        SQUE_NR = ncolors;
      n = ncolors;
      ak = 105.; s = 150.; c = s/60.;
      ab = 65.; sb = 190.; cb = s/60.;
      m=256;
      an = a1;
      da = a2 - a1;     da = da/n;
      an = an-da+360.;

      for(i=0; i < ncolors; i++) {
        an += da;   an=fmod(an,360.);
        if((an >= 0) && (an < 120.)) {
          r = 255.*pow((ak + min(s,(120. - an)*c))/255., gamm) +.5;
          g = 255.*pow((ak + min(s,an*c))/255., gamm) +.5;
          b=0;
        }
        if((an >= 120.) && (an < 240.)) {
          r = 0;
          g = 255.*pow((ak + min(s ,(240. - an)*c))/255., gamm) +.5;
          b = 255.*pow((ab + min(sb,(an - 120.)*cb))/255., gamm) +.5;
        }
        if(an >= 240.) {
          r =  255.*pow((ak + min(s,(an - 240.)*c))/255., gamm) +.5;
          g = 0;
          b = 255.*pow((ak + min(s,(360. - an)*c))/255., gamm) +.5;
        }
        MYcol[i].pixel = pixels[ncolors-1-i];
        MYcol[i].red   = r*m;
        MYcol[i].green = g*m;
        MYcol[i].blue  = b*m;
        MYcol[i].flags = DoRed|DoGreen|DoBlue;
      }
      No_Color_init = 0;
   }
   XStoreColors(theDisp, CMap, MYcol, ncolors);
}

/* --------------- */
   c_rotate(k)
   int k;
/* --------------- */
{
   if (YES_color) color_rotate(k);
   else           grey_rotate(-k);
}


/* --------------- */
   color_rotate(k)
   int k;
/* --------------- */
{
   register int i, j;

   FIRST = 1;
   SQUE_NR = ncolors;


   if(k >= 0) {
     for(i=0; i < k; i++) {
     tmp1[i] = MYcol[i].red;
     tmp2[i] = MYcol[i].green;
     tmp3[i] = MYcol[i].blue;
     }
     for(i=0; i < ncolors - k; i++) {
       MYcol[i].red   = MYcol[i+k].red;
       MYcol[i].green = MYcol[i+k].green;
       MYcol[i].blue  = MYcol[i+k].blue;
     }
     for(i = ncolors - k, j = 0; i < ncolors; i++, j++) {
       MYcol[i].red   = tmp1[j];
       MYcol[i].green = tmp2[j];
       MYcol[i].blue  = tmp3[j];
     }
   }
   else {
     k = -k;
     for(i = ncolors - 1; i >= ncolors - k; i--) {
     tmp1[i] = MYcol[i].red;
     tmp2[i] = MYcol[i].green;
     tmp3[i] = MYcol[i].blue;
     }
     for(i = ncolors - 1; i >= k; i--) {
       MYcol[i].red   = MYcol[i-k].red;
       MYcol[i].green = MYcol[i-k].green;
       MYcol[i].blue  = MYcol[i-k].blue;
     }
     for(i = k-1, j = ncolors - 1; i >=0; i--, j--) {
       MYcol[i].red   = tmp1[j];
       MYcol[i].green = tmp2[j];
       MYcol[i].blue  = tmp3[j];
     }
   }
   XStoreColors(theDisp, CMap, MYcol, ncolors);
}

/* --------------- */
   grey_rotate(k)
   int k;
/* --------------- */
{
   register int i, j;

   FIRST = 1;
   SQUE_NR = ncolors;


   if(k >= 0) {
     for(i=0; i < k; i++) {
     tmp1[i] = MYgrey[i].red;
     tmp2[i] = MYgrey[i].green;
     tmp3[i] = MYgrey[i].blue;
     }
     for(i=0; i < ncolors - k; i++) {
       MYgrey[i].red   = MYgrey[i+k].red;
       MYgrey[i].green = MYgrey[i+k].green;
       MYgrey[i].blue  = MYgrey[i+k].blue;
     }
     for(i = ncolors - k, j = 0; i < ncolors; i++, j++) {
       MYgrey[i].red   = tmp1[j];
       MYgrey[i].green = tmp2[j];
       MYgrey[i].blue  = tmp3[j];
     }
   }
   else {
     k = -k;
     for(i = ncolors - 1; i >= ncolors - k; i--) {
     tmp1[i] = MYgrey[i].red;
     tmp2[i] = MYgrey[i].green;
     tmp3[i] = MYgrey[i].blue;
     }
     for(i = ncolors - 1; i >= k; i--) {
       MYgrey[i].red   = MYgrey[i-k].red;
       MYgrey[i].green = MYgrey[i-k].green;
       MYgrey[i].blue  = MYgrey[i-k].blue;
     }
     for(i = k-1, j = ncolors - 1; i >=0; i--, j--) {
       MYgrey[i].red   = tmp1[j];
       MYgrey[i].green = tmp2[j];
       MYgrey[i].blue  = tmp3[j];
     }
   }
   XStoreColors(theDisp, CMap, MYgrey, ncolors);
}

/* -------- */
   c_swap()
/* -------- */
{
   if (YES_color) color_swap();
}

/* --------------- */
   color_swap()
/* --------------- */
{
   register int i, k;

   FIRST = 1;	SQUE_NR = ncolors;
   k = ncolors - 1;

   for(i=0; i < ncolors; i++) {
    tmp1[i] = MYcol[i].red;
    tmp2[i] = MYcol[i].green;
    tmp3[i] = MYcol[i].blue;
   }
   for(i=0; i < ncolors ; i++) {
     MYcol[i].red   = tmp1[k-i];
     MYcol[i].green = tmp2[k-i];
     MYcol[i].blue  = tmp3[k-i];
   }
   XStoreColors(theDisp, CMap, MYcol, ncolors);
}

/* ------------- */	/* Modify span or contrast */
   c_squeeze(d)
   int d;
/* ------------- */
{
   if (YES_color) color_squeeze(d);
   else           grey_contrast(-2*d);
}

/* -------------------- */      /* Modify contrast of the Image */
   grey_contrast(d_lev)
   int d_lev;
/* -------------------- */
{
   register int i, k, m, delta, dx;

   dx = max(((abs(INgrey[ncolors-1] - INgrey[0])>>6)/ncolors), 1);
   delta = d_lev*dx;
   for(i=0; i < ncolors; i++) {
      k = INgrey[i] += i*delta;
      m = min_max_col(k);
      MYgrey[i].red   = m;
      MYgrey[i].green = m;
      MYgrey[i].blue  = m;
   }
   XStoreColors(theDisp, CMap, MYgrey, ncolors);       /* Set for printer */
}

/* -------------------- */	/* Modify span of colors*/
   color_squeeze(d)
   int d;
/* -------------------- */
{
   static unsigned int x_arr[MCOLORS];
   register int i, ex, w;

   w = ncolors/2;
   SQUE_NR = max(1, SQUE_NR - d);
   for (ex = 0; ex < ncolors; ex++)  x_arr[ex] = SQUE_NR*(ex - w)/ncolors + w;

   if(FIRST)
   for(i=0; i < ncolors; i++) {
     tmp1[i] = MYcol[i].red;
     tmp2[i] = MYcol[i].green;
     tmp3[i] = MYcol[i].blue;
     FIRST = 0;
   }
   for(i=0; i < ncolors; i++) {
     if(x_arr[i] > uucolors) {
       MYcol[i].red   = 0;
       MYcol[i].green = 0;
       MYcol[i].blue  = 0;
     }
     else {
       MYcol[i].red   = tmp1[x_arr[i]];
       MYcol[i].green = tmp2[x_arr[i]];
       MYcol[i].blue  = tmp3[x_arr[i]];
    }
  }
  XStoreColors(theDisp, CMap, MYcol, ncolors);
}

/* ------------------ */	/* Modify brightness (color or B&W) */
   c_bright(d)
   int d;
/* ------------------ */
{
   if (YES_color) color_bright(d);
   else           grey_change(-2*d);
}

/* ------------------ */   /* Modify brightness of the Image  */
   grey_change(d_lev)
   int d_lev;
/* ------------------ */
{
   register int i, k, m, delta, dx;

   dx = abs((INgrey[ncolors-1] - INgrey[0])/ncolors);
   delta = d_lev*dx;
   for (i=0; i < ncolors; i++) {
      k = INgrey[i] += delta;
      m = min_max_col(k);
      MYgrey[i].red   = m;
      MYgrey[i].green = m;
      MYgrey[i].blue  = m;
   }
   XStoreColors(theDisp, CMap, MYgrey, ncolors);       /* Set for printer */
}

/* ------------------ */	/* Modify brightness of the Image	*/
   color_bright(d)
   int d;
/* ------------------ */
{
   double c;
   register int i;
   c = d; c = 1. - c*.005;
   FIRST = 1;	SQUE_NR = ncolors;
	
   for(i=0; i < ncolors; i++) {
     MYcol[i].red   = min_max_col(c*MYcol[i].red);
     MYcol[i].green = min_max_col(c*MYcol[i].green);
     MYcol[i].blue  = min_max_col(c*MYcol[i].blue);
   }
   XStoreColors(theDisp, CMap, MYcol, ncolors);
}

/* -------------- */ /* It assigns index for one of standard 16 colors */
   STD_colors(cx)    /* First index is 1, last 16 (equv. to 0 - black col.). */
   short int cx;
/* -------------- */
{
   XColor  any_col;

   if ( color_x11[cx-1] + 1 ) return(color_x11[cx-1]);
   else {
      any_col.red   = Solid_color[cx-1].red;
      any_col.green = Solid_color[cx-1].green;
      any_col.blue  = Solid_color[cx-1].blue;
      any_col.flags = DoRed | DoGreen | DoBlue;
      if ( !XAllocColor(theDisp, CMap, &any_col) )
         FatalError ("XAllocColor problem in STD_colors(). ");
      return( color_x11[cx-1] = any_col.pixel );
   }
}

/* ---------------------------------- */   /* Create Image from the im_arr */
   XImage *Load_Any_Arr(im_arr, x, y)      /* Usage: theImage = Load_An... */
   short int im_arr[];                     /* Uses own 16 colors for nega- */
   int       x, y;                         /* tive numbers (-1 - -16).     */
/* ---------------------------------- */   /* 8 and 12 bit planes work OK. */
{
  register char *ptr;
  register int   i, j, k, iN, iE, w2;
  int        Width, Hight;
  char      *Image;
  XImage    *image;

  iN = (Planes + 7)/8;             /* 1 or 2, pointer increment */
  iE = iN - 1;                     /* 0 or 1 for next pointer */

  Width = x;      /* Image width  */
  Hight = y;      /* Image higth  */
  w2 = Width * iN;

  Image = (char *) malloc(w2*Hight);
  if (!Image) FatalError("Not enough memory for the Image");
  image = XCreateImage(theDisp, theVisual, Planes, ZPixmap, 0, Image,
                                Width, Hight, 8, w2);
  if (! image) FatalError("Unable create image in Load_Any_Arr().");

  ptr = Image;
  k = 0;
  for (i=0; i<Hight; i++)
  for (j=0; j<Width; j++, ptr += iN)
  if(im_arr[k] >= 0) {
    *ptr = pixels[im_arr[k]] >> 8;
    *(ptr + iE) = pixels[im_arr[k++]];
  }
  else {
    *ptr = STD_colors(-im_arr[k]) >> 8;
    *(ptr + iE) = STD_colors(-im_arr[k++]);
  }
  return image;
}

/* -------------------------------------------- */
   New_Cursor(Disp, Wind, shape, fg_col, bg_col)
   Display	*Disp;
   Window	Wind;
   unsigned int	shape;
   char 	*fg_col, *bg_col;
/* -------------------------------------------- */
/* Assigns new cursor to the Window. Call after creating the Window */
{
   XColor	Im_cur_fg, Im_cur_bg;
   Colormap	CM;
   Cursor	cursor;

   CM = DefaultColormap(Disp, DefaultScreen(Disp));

   cursor = XCreateFontCursor(Disp, shape);        /* new cursor shape */
   XDefineCursor(Disp, Wind, cursor);
   XParseColor(Disp, CM, fg_col, &Im_cur_fg);    /* forground  color */ 
   XParseColor(Disp, CM, bg_col, &Im_cur_bg);    /* background color */
   XRecolorCursor(Disp, cursor, &Im_cur_fg, &Im_cur_bg);
}

/* ----------------------------------------------------------------
   Read file subroutine fo use in C.        A.Jesmanowicz, MCW 1991
	return error :	0 - OK,
			1 - opening problem,
			2 - file longer then array.
	fname :	file name.
	size  :	on input - max size of the arr or 0 for any length,
		on output- real size of the file (and arr in bytes).
	arr   :	returned file as array.
   ---------------------------------------------------------------- */
/* ----------------------------- */
   int  read_iqm(fname,size,arr)
   int  *size;
   char fname[],arr[];
/* ----------------------------- */
{
	int	isize = *size;
	int	fp;				/* file descriptor  */
	struct	stat file_stat;			/* status structure */

	if ((fp = open(fname, O_RDONLY)) <= 0)	/* file must exist */
           return(1); 		             	/* or error = 1.   */

	fstat(fp, &file_stat);			/* get file size in bytes   */

	if(file_stat.st_size > isize && isize)	/* file can not be too long */
	   return(2);     	                /* or error = 2.	    */

	*size =  file_stat.st_size;		/* return file size */

	read(fp, arr, file_stat.st_size);	/* read whole file  */
	close(fp);
	return(0);                              /* no error : 0 */
}

/* ----------------------------------------------------------------
   Write file subroutine fo use in C.       A.Jesmanowicz, MCW 1991
	return error :	0 - OK,
			1 - opening problem,
			2 - negative length.
	fname :	file name.
	size  :	size of the file.
	arr   :	file array.
   ---------------------------------------------------------------- */
/* ----------------------------- */
   int WRite_iqm(fname,size,arr)
   int  *size;
   char fname[],arr[];
/* ----------------------------- */
{
	int	isize = *size;
	int	fp;				/* file descriptor  */

	if(isize < 0)				/* size has to be real */
	   return(2);				/* or error = 2.       */

	if ((fp = open(fname, O_WRONLY|O_CREAT|O_EXCL,0644)) <= 0)
	   return(1);				/* file must not exist */
						/* or error = 1.       */

	write(fp, arr, isize);			/* write whole file */
	close(fp);
	return(0);
}


/*  Directory names (if present) are skipped. and number is added in front */
/* -------------- */
   top_nr_name(name)
   char *name;
/* -------------- */
{
   register int i, k;
   char         *p_name;

   k = strlen( name );
   p_name = name;
   for(i = 0; i < k; i++) {
      if(name[i] != 32 && name[i]) {
         if(name[i] == 47) p_name = &name[i+1];
      }
      else {
         name[i] = 0;
         break;
      }
   }
   sprintf(T_name, "#%d   %s", Im_Nr+1, p_name);
}

/* ------------------ */ /* resample image into 256x256 frame using st_? */
   void Resample(all)
   short int *all;
/* ------------------ */
{
   register short int *a;
   register int       *s, i1, ix, lx, i, j, k, l, kk, mm, dkk, ii, ij, m, n;
   register int       II, J, IJ, IJk, ijl, k0, IM2, lx2, ly, IJK;

      a = all;

      for (i=0; i< 256*256; i++)
         tmp_imx[i] = imx[i] = 0;     /* zero everything to begin with AWS */

      for (i=0; i < 256*16; i++)
       {  m = i%64 + 256*(63 - i/64);
         tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*16; i < 256*32; i++)
       {  m = i%64 + 64 + 256*(63-(i-256*16)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*32; i < 256*48; i++)
       {  m = i%64 + 128 + 256*(63-(i-256*32)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*48; i < 256*64; i++)
       {  m = i%64 + 192 + 256*(63-(i-256*48)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*64; i < 256*80; i++)
       {  m = i%64 + 256*64 + 256*(63-(i-256*64)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*80; i < 256*96; i++)
       {  m = i%64 + 64 + 256*64 + 256*(63-(i-256*80)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*96; i < 256*112; i++)
       {  m = i%64 + 128 +256*64 + 256*(63-(i-256*96)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*112; i < 256*128; i++)
       {  m = i%64 + 192 + 256*64 + 256*(63-(i-256*112)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*128; i < 256*144; i++)
       {  m = i%64 + 256*128 + 256*(63-(i-256*128)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*144; i < 256*160; i++)
       {  m = i%64 + 64 + 256*128 + 256*(63-(i-256*144)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*160; i < 256*176; i++)
       {  m = i%64 + 128 + 256*128 + 256*(63-(i-256*160)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*176; i < 256*192; i++)
       {  m = i%64 + 192 + 256*128 + 256*(63-(i-256*176)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*192; i < 256*208; i++)
       {  m = i%64 + 256*192 + 256*(63-(i-256*192)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*208; i < 256*224; i++)
       {  m = i%64 + 64 + 256*192 + 256*(63-(i-256*208)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*224; i < 256*240; i++)
       {  m = i%64 + 128 + 256*192 + 256*(63-(i-256*224)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      for (i=256*240; i < 256*256; i++)
       {  m = i%64 + 192 + 256*192 + 256*(63-(i-256*240)/64);
          tmp_imx[m] = imx[m] =  a[i]; }

      return;
}

/* ------------ */ /* Don't use before first Resample() and Load_Any_Arr(). */
   void            /* Load_Any_Arr() allocate memory for theImage. */
   Put_image(n)    /* Put_image() needs THIS memory. */
   int n;
/* ------------ */
{
   register int    min2, max2, i, m;
   register float  coef2, del2;
   int             xtemp1, xtemp2, ytemp1, ytemp2, xmat, xm, ym, xxx;

   if ( n == old_im ) { 
      for (i = 0; i < IM_ARR; i++) imx[i] = tmp_imx[i];
   }
   else {
      m = n * ar_size;
      if ( n > (npoints - 1) ) {   /* Extra images, out of time course */ 
         min2 =  32767;
         max2 = -32768;
         for (i=0; i < ar_size; i++) {
            if ( ar[m+i] < min2 ) min2 = ar[m+i];
            if ( ar[m+i] > max2 ) max2 = ar[m+i];
         }
         del2 = max2 - min2;
         if (del2 < 1.) del2 = 1.;
         coef2 = ( (float) NC - 1.) / del2;
      }
      else {
         coef2 = coef1;
         min2 = min1;
      }
                             /* make differential single image */
      if ( diff_im || avr_grp ) {
         if ( ! avr_grp )
            for (i=0; i < ar_size; i++) tmp_ar[i+offs] = ar[m+i] - ref_ar[i];

                              /* make average image */
         else {
            Av_1 = n;
            Av_2 = Av_1 + Av_length - 1;
            for (i=0; i< ar_size; i++) av_ar[i] = 0;
            for (i=0; i < ar_size; i++) {
               for (j=Av_1; j <= Av_2; j++)
                  av_ar[i] += ar[j*ar_size+i];
            }
            for (i=0; i < ar_size; i++)  av_ar[i] /= Av_length;

                              /* average against ref image */
            if ( avr_grp && diff_im )
               for (i=0; i < ar_size; i++)
                  tmp_ar[i+offs] = av_ar[i] - ref_ar[i];

                              /* pure average image */
            else 
               for (i=0; i < ar_size; i++)
                  tmp_ar[i+offs] = av_ar[i];

         }
                              /* renormalize image */
         min2 =  32767;
         max2 = -32768;
         for (i=0; i < ar_size; i++) {
            if ( tmp_ar[i+offs] < min2 ) min2 = tmp_ar[i+offs];
            if ( tmp_ar[i+offs] > max2 ) max2 = tmp_ar[i+offs];
         }
         del2 = max2 - min2;
         if (del2 < 1.)  del2 = 1.;
         coef2 = ( (float) NC - 1.) / del2;
         for (i=0; i < ar_size; i++)
            tmp_ar[i+offs] = (int)((float)(tmp_ar[i+offs]-min2)*coef2 + .5);
      }
      else                                        /* Normal Image */
         for (i=0; i < ar_size; i++)              /* renormalize data */
            tmp_ar[i+offs] = (int)((float)(ar[m+i]-min2)*coef2 + .5);
   
      Resample(&tmp_ar[offs]);
      top_nr_name(f_name[n]);
      Allow_new_name(Argc, Argv);
      DrawTopWindow();
      old_im = n;
   }

   xxx = max(0, x_mag/2 - 1);        /* stuff for proper pivot and frame */
   xm = x_mag - 1;
   ym = x_mag - 1;
   xmat = x_mag * (mat - 1) + x_mag/2 + 1;
   
   ytemp1 = x_mag * (ypoint + yc) + ym ;
   if (ytemp1 < 0)          ytemp1 += IM_HEIGHT;
   if (ytemp1 >= IM_HEIGHT) ytemp1 -= IM_HEIGHT;

   ytemp2 = x_mag * (ypoint + yc - mat + 1);
   if (ytemp2 < 0)          ytemp2 += IM_HEIGHT;
   if (ytemp2 >= IM_HEIGHT) ytemp2 -= IM_HEIGHT;

   for (i=0; i < xmat + xxx; i++) {
     xtemp1 = x_mag * (xpoint - xc) + i;
     if (xtemp1 < 0)          xtemp1 += IM_HEIGHT;
     if (xtemp1 >= IM_HEIGHT) xtemp1 -= IM_HEIGHT;

     imx[ytemp1 * IM_HEIGHT + xtemp1] = fr_index;
     imx[ytemp2 * IM_HEIGHT + xtemp1] = fr_index;
   }

   xtemp1 = x_mag * (xpoint - xc);
   if (xtemp1 < 0)          xtemp1 += IM_HEIGHT;
   if (xtemp1 >= IM_HEIGHT) xtemp1 -= IM_HEIGHT;

   xtemp2 = x_mag * (xpoint - xc + mat - 1) + xm;
   if (xtemp2 < 0)          xtemp2 += IM_HEIGHT;
   if (xtemp2 >= IM_HEIGHT) xtemp2 -= IM_HEIGHT;

   for (i=0; i < xmat; i++) {
      ytemp1 = x_mag * (ypoint + yc) - i + ym;
      if (ytemp1 < 0)          ytemp1 += IM_HEIGHT;
      if (ytemp1 >= IM_HEIGHT) ytemp1 -= IM_HEIGHT;
   
      imx[ytemp1 * IM_HEIGHT + xtemp1] = fr_index;
      imx[ytemp1 * IM_HEIGHT + xtemp2] = fr_index;
   }

   Load_Next_Arr(theImage, imx,IM_HEIGHT,IM_HEIGHT);
   MResize(&expImage,&eWIDE,&eHIGH,theImage,eWIDE,eHIGH);
   XPutImage(theDisp, theWindow, theGC, expImage,0, 0, 0, 0, eWIDE, eHIGH);
}

/* ---------------------------------- */   /* Create Image from the im_arr */
   Load_Next_Arr(Image, im_arr, x, y)      /* Usage: Load_Next... */
   short int im_arr[];                     /* Uses own 16 colors for nega- */
   int       x, y;                         /* tive numbers (-1 - -16).     */
   XImage    *Image;
/* ---------------------------------- */   /* 8 and 12 bit planes work OK. */
{
  register char *ptr;
  register int   i, j, k, iN, iE;
  int        Width, Hight;

  iN = (Planes + 7)/8;             /* 1 or 2, pointer increment */
  iE = iN - 1;                     /* 0 or 1 for next pointer */

  Width = x;      /* Image width  */
  Hight = y;      /* Image higth  */

  ptr = Image->data;
  k = 0;
  for (i=0; i<Hight; i++)
  for (j=0; j<Width; j++, ptr += iN)
  if(im_arr[k] >= 0) {
    *ptr = pixels[im_arr[k]] >> 8;
    *(ptr + iE) = pixels[im_arr[k++]];
  }
  else {
    *ptr = STD_colors(-im_arr[k]) >> 8;
    *(ptr + iE) = STD_colors(-im_arr[k++]);
  }
}
/* ---------------------------- */
   Allow_new_name(argc, argv)
   int  argc;
   char *argv[];
/* ---------------------------- */
{
   XSetStandardProperties(theDisp, theWindow, T_name, T_name, None,
                          argv, argc, NULL);
}

/* ----------------------------- */   /* plot scaled new_line to pixmap  */
   void plot_line()
/* ----------------------------- */
{
   register int i, m, temp, index, ix, iy, xtemp, ytemp;

   m = IM_ARR;

   for (ix=0;ix<mat;ix++) {
      xtemp = xpoint + ix - xc;
      if (xtemp < 0)        xtemp += im_size;
      if (xtemp >= im_size) xtemp -= im_size;
      for (iy=0;iy<mat;iy++) {
         ytemp = ypoint - iy + yc;
         if (ytemp < 0)        ytemp += im_size;
         if (ytemp >= im_size) ytemp -= im_size;
         index = ytemp * im_size + xtemp;
         for (i=0; i < N_im; i++) val[ix][iy][i] = zzz[i*m+index];
         pmin[ix][iy] = 40000;
         for (i=0; i < npoints; i++)
            pmin[ix][iy] = min(pmin[ix][iy],val[ix][iy][i]);
         if (iscale > 0) {
           for (i=0; i < npoints; i++)
              plot[ix][iy][i] = (val[ix][iy][i]-pmin[ix][iy])*iscale;
           pmax[ix][iy] = pmin[ix][iy] + gy / iscale;
         }
         else {
            temp = -iscale;
            for (i=0; i < npoints; i++)
               plot[ix][iy][i] = (val[ix][iy][i]-pmin[ix][iy])/temp;
            pmax[ix][iy] = pmin[ix][iy] + gy * temp;
         }
         line_color("black");

         for (i=0; i < npoints; i++) {
            a_line[i].x = xorigin[ix][iy] + i*gx/(npoints-1);
            a_line[i].y = iHIGH - yorigin[ix][iy] - plot[ix][iy][i];
         XFillRectangle(theDisp, pxWind, theGC, a_line[i].x, a_line[i].y, jimpoint, jimpoint);
         }

         XDrawLines(theDisp, pxWind, theGC, a_line, npoints, CoordModeOrigin);
      }
   }
}

/* ----------------------------- */ /* draw marker for current image  */
   void draw_marker()
/* ----------------------------- */ /* marker fixed (shifted by one im)  */
{
   register int i, j, k, g, xo, yo, x1, dx;

   if (mark) {
     line_color(color[color_index]);
     g = grid_ar[grid_index];
     for (i=0;i<mat;i++) {
        for (j=1;j<=npoints/g;j++) {
	      k = xorigin[i][0] + (j * g - 1) * gx / (npoints - 1);
          plotx(k,mdy1, 0);
          plotx(k,mdy1+mat*gy, 1);
	    }
     }
   }

   line_color("white");
   xo = xorigin[xc][yc];
   yo = iHIGH - yorigin[xc][yc] - gy;
   XFillRectangle(theDisp, pxWind, theGC, xo, yo, gx, gy);

   if (diff_im) {
      line_color("blue");
      x1 = Im_1*gx/(npoints-1);
      dx = (Im_2 - Im_1)*gx/(npoints-1) + 1;
      XFillRectangle(theDisp, pxWind, theGC, xo + x1, yo, dx, gy);
   }

 /*  if (avr_grp) {
      line_color("red");
      x1 = Av_1*gx/(npoints-1);
      dx = (Av_2 - Av_1)*gx/(npoints-1) + 1;
      XFillRectangle(theDisp, pxWind, theGC, xo + x1, yo, dx, gy);
   }  */

   if (mark) {
     line_color("white");
     for (j=1;j<=npoints/g;j++) {
	   k = xorigin[xc][yc] + (j * g - 1) * gx / (npoints - 1);
       plotx(k,yorigin[xc][yc]   , 0);
       plotx(k,yorigin[xc][yc]+gy, 1);
     }
   }
}

/* ----------------------------- */   /* scale plot up and redraw  */
   void scale_up()
/* ----------------------------- */
{
   if (iscale > 0) iscale *= 2;
   else if (iscale < -2) iscale /= 2;
   else iscale = 1;
}

/* ----------------------------- */   /* scale plot up and redraw  */
   void scale_down()
/* ----------------------------- */
{
   if (iscale > 1) iscale /= 2;
   else if (iscale < 0) iscale *= 2;
   else iscale = -2;
}

/* ----------------------------- */   /* decrease matrix and redraw  */
   void mat_down()
/* ----------------------------- */
{
   int old;
   
   old = mat;
   mat--;
   if (mat < 1) mat = 1;
   else if (mat > MAT_MAX) mat = MAT_MAX;
   if (mat!= old) {
      init_mat();
      Put_image(Im_Nr);
      redraw_graph();
      DrawSubWindow();
   }
}

/* ----------------------------- */   /* increase matrix and redraw  */
   void mat_up()
/* ----------------------------- */
{
   int old;
   
   old = mat;
   mat++;
   if (mat < 1) mat = 1;
   else if (mat > MAT_MAX) mat = MAT_MAX;
   if (mat!= old) {
      init_mat();
      Put_image(Im_Nr);
      redraw_graph();
      DrawSubWindow();
   }
}

/* ----------------------------- */   /* initialize matrix stuff  */
   void init_mat()
/* ----------------------------- */
{
   gx = GX_MAX / mat;
   gy = GY_MAX / mat;
   for (i=0;i<mat;i++) {
     for (j=0;j<mat;j++) {
	   xorigin[i][j] = mdx1 + i * gx;
	   yorigin[i][j] = mdy1 + j * gy;
	 }
   }
   xc = mat/2;
   yc = (mat-1)/2;
}

/* ----------------------------- */   /* decrease grid spacing and redraw  */
   void grid_down()
/* ----------------------------- */
{
   int old;
   
   old = grid_index;   
   grid_index--;
   if (grid_index < 0) grid_index = 0;
   if (grid_index!= old) {
      redraw_graph();
      DrawSubWindow();
   }
}

/* ----------------------------- */   /* increase grid spacing and redraw  */
   void grid_up()
/* ----------------------------- */
{
   int old;
   
   old = grid_index;   
   grid_index++;
   if (grid_index >= GRID_NUM) grid_index = GRID_NUM - 1;
   if (grid_index!= old) {
      redraw_graph();
      DrawSubWindow();
   }
}

/* ------------ */
   init_const()
/* ------------ */
{
   iscale = 4;
   auto_scale = .4;

   xpoint = im_size/2;
   ypoint = im_size/2;

   xspace = 5;
   yspace = 20;
   mytxt = 20;
   mdx1   = GL_DLX + 1;
   mdy1   = GB_DLY + 1;
   idx   = GL_DLX + GX_MAX + GR_DLX;
   idy   = GB_DLY + GY_MAX + GT_DLY;
   grid_ar[0] = 2;
   grid_ar[1] = 5;
   grid_ar[2] = 10;
   grid_ar[3] = 20;
   grid_ar[4] = 50;
   grid_index = 2;
   mark = 1;
   color_index = 0;
   mat = 5;
   init_mat();
}

/* ------------------------------------------------- Main link from main() C */
   void window_plane()
/* --------------------------------------------------- */
{
   int  i, out_loop = 0, tmpx;

   XGCValues   gcv;

   iWIDE = idx;	iHIGH = idy;

   CreateGraphWindow(NULL, 0);      /* Set theWindow properties */
/* subWindow, topWindow and keys stuff is here ------- vvvvvvv ------- AJ*/

    kfont = XGetDefault(theDisp, ProgramName, "KeyFont");
    if (!kfont) kfont = KFONT;

    ffc = XGetDefault(theDisp, ProgramName, "FKeyFore");
    fbc = XGetDefault(theDisp, ProgramName, "FKeyBack");

    /* Set normal default colors */
    ForeColor = BlackPixel(theDisp, theScreen);
    BackColor = WhitePixel(theDisp, theScreen);

    FKeyFore =  ForeColor;
    FKeyBack =  BackColor;

    if ( (ffc != NULL) && XAllocNamedColor(theDisp, CMap, ffc, &xcsd, &xced));
    else if (XAllocNamedColor(theDisp, CMap, "white", &xcsd, &xced));
    else    FatalError ("XAllocNamedColor problem. ");
    FKeyFore=xcsd.pixel;


    if ( (fbc != NULL) && XAllocNamedColor(theDisp, CMap, fbc, &xcsd, &xced));
    else if (XAllocNamedColor(theDisp, CMap, "blue", &xcsd, &xced));
    else    FatalError ("XAllocNamedColor problem. ");
    FKeyBack=xcsd.pixel;

    /* load fonts, figure out sizes of keypad and display */
    if ( !(kfontinfo = XLoadQueryFont(theDisp, kfont)) )
       if ( !(kfontinfo = XLoadQueryFont(theDisp, "fixed")) ) {
          sprintf(strp, "Can't open %s or fixed fonts\n", kfont);
          FatalError (strp);
       }
    keyfont = kfontinfo->fid;

    minwide = XTextWidth(kfontinfo, "MMMM", 4);
    keyhigh = kfontinfo->ascent + 7;
    for (i=0; i<N_KEYS; i++) {
       tmpx = XTextWidth(kfontinfo, xtkeys[i].st, strlen(xtkeys[i].st));
       keywide[i] = max(minwide, tmpx) + 2;
    }

    gcv.foreground = FKeyFore;
    gcv.function = GXcopy;
    gcv.font = keyfont;
    Fkeygc = XCreateGC(theDisp, GWindow, GCForeground|GCFunction|GCFont, &gcv);

    gcv.function = GXinvert;
    gcv.plane_mask = FKeyFore ^ FKeyBack;
    Fkeyigc = XCreateGC(theDisp, GWindow, GCFunction|GCPlaneMask, &gcv);

    Setup_subWindow();

    Setup_topWindow();

    Setup_keys();

/* subWindow, topWindow and keys stuff is here ------- ^^^^^^^ ------- AJ*/

                         /* Make dummy window for text and graphic entry */
   pxWind = XCreatePixmap(theDisp, GWindow, iWIDE, iHIGH, Planes);

   New_Cursor(theDisp, GWindow, XC_left_ptr, "blue", "yellow"); /* cursor */

   /* Text color and font */

   XSelectInput(theDisp, GWindow, ExposureMask | KeyPressMask
                 | ButtonPressMask);
   XMapWindow(theDisp,GWindow);		/* Show window first time */
   XMapSubwindows(theDisp, GWindow);    /* All keys stuff */

   
   for (i = 0; i < N_KEYS+2; i++) exp_done[i] = 0;

   while (1) {
      Window wind;

      XNextEvent(theDisp, &first_event);
                                          /* Wait for window on */
      switch (first_event.type) {
         case Expose: {
            XExposeEvent *exp_event = (XExposeEvent *) &first_event;
            wind = exp_event->window;

            if (wind ==   GWindow)
               exp_done[N_KEYS] = 1;
            else if (wind == subWindow) {
               DrawSubWindow();
               exp_done[N_KEYS+1] = 1;
            }
            else if (wind == topWindow) {
               DrawTopWindow();
               exp_done[N_KEYS+2] = 1;
            }
            else {
               for (i=0; i < N_KEYS; i++) {
                  if (wind == key[i].wid) {
                     if (i > 2) DrawKey(i);
                     exp_done[i] = 1;
                  }
               }
            }
            out_loop = exp_done[N_KEYS] * exp_done[N_KEYS+1]
                                        * exp_done[N_KEYS+2];
            for (i = 0; i < N_KEYS; i++)  out_loop *= exp_done[i];
            if (out_loop) return;
         }
         break;

         default:                /* ignore unexpected events */
         break;
      }
   }
}

/* ----------------- */
   Setup_subWindow()
/* ----------------- */
{
   unsigned long white;

   if (!(XAllocNamedColor(theDisp, CMap, "white", &any_col, &rgb_col)))
   FatalError ("XAllocNamedColor problem. AJ in Setup_subWindow()");
   white = any_col.pixel;

   sub_W_x = idx;
   sub_W_y = GB_DLY;
   subWindow = XCreateSimpleWindow(theDisp, GWindow, 0, GY_MAX + GT_DLY,
                     sub_W_x, sub_W_y, 1, white, white);
   XSelectInput(theDisp, subWindow, ExposureMask | ButtonPressMask |
                  ButtonReleaseMask);
}

/* --------------- */ /* Draw subWindow containing text info */
   void
   DrawSubWindow()
/* --------------- */
{
   float   diff_prcnt, stdev;
   int     str_x = 5, x, y, loc, i, kkk, kkkx, kkky;
   int     zzzz1, temph, templ, imnh, imnl, mean, mean2, ptop;

   line_color("blue");
   XFillRectangle(theDisp, subWindow, theGC, 0, 0, idx, GB_DLY);

/*   sprintf (strp,"Im. rot.");
   subW_TXT(str_x, 35, strp);
   sprintf (strp,"clockwise");
   subW_TXT(str_x, 20, strp);
   sprintf(strp, "%d*90 deg.", rot_nr);
   subW_TXT(str_x,  5, strp);
   str_x += XTextWidth(mfinfo ,"clockwise", strlen("clockwise")) + 30;  */

   txt_color("yellow");

   kkkx = xpoint % 64 + 1;
   kkky = 64 - ypoint % 64;
   kkk = 4 * (ypoint/64) + xpoint/64 + 1;
   sprintf (strp,"X: %d", kkkx);
   subW_TXT(GL_DLX + 10, 28, strp);
   sprintf (strp,"Y: %d", kkky);
   subW_TXT(GL_DLX + 10, 13, strp);
   sprintf (strp, "AT  SLICE");
   subW_TXT(GL_DLX + 50, 20, strp);
   sprintf (strp, "%d", kkk);
   subW_TXT(GL_DLX + 120, 20, strp);

   str_x += XTextWidth(mfinfo, "X: 000", strlen("X: 000")) + 170;

   sprintf (strp,"Pix. value: %d at vol. #%d", val[xc][yc][Im_Nr], Im_Nr+1);
   subW_TXT(str_x, 20, strp);
   loc = str_x;

   temph = 0;
   templ = 32767;
   imnh = 0;
   imnl = 0;
   mean = 0;
   mean2 = 0;
   stdev = 0.0;

   for (zzzz1 = 0; zzzz1 < npoints; zzzz1++)
   {    if (temph < val[xc][yc][zzzz1]) {
              temph = val[xc][yc][zzzz1];
              imnh  = zzzz1 + 1;  }
        if (templ > val[xc][yc][zzzz1]) {
              templ = val[xc][yc][zzzz1];
              imnl = zzzz1 + 1; }
        mean += val[xc][yc][zzzz1];
        mean2 += val[xc][yc][zzzz1] * val[xc][yc][zzzz1];
   }

   ptop = temph -templ;
   stdev = (float)(sqrt((npoints*mean2 - mean* mean)/(npoints*(npoints-1))));
   mean = mean / npoints;


   sprintf (strp,"peak max %d at vol. #%d,  peak min %d at vol. #%d", temph, imnh, templ, imnl);
   subW_TXT(GL_DLX + 330, 28, strp);

   sprintf (strp,"mean %d,  standard deviation %3.1f,  peak2peak %d", mean, stdev, ptop);
   subW_TXT(GL_DLX + 330, 13, strp);


 /*  sprintf (strp,"Grid:%d", grid_ar[grid_index]);
   subW_TXT(GL_DLX + 360, 20, strp);   */

   sprintf(strp, "allen@iris.nimh.nih.gov");
   subW_TXT(GL_DLX + 640, 5, strp);

   sprintf(strp, "Total points: %d", npoints);
   subW_TXT(GL_DLX + 640, 20, strp);
   txt_color("blue");

   i = ypoint*im_size + xpoint;
   if ( diff_im) {
      if ( avr_grp ) {
         diff_prcnt = (float) (200 * (av_ar[i] - ref_ar[i]))
                           / (float) (av_ar[i] + ref_ar[i]);
      }
      else {
         diff_prcnt = (float) (200 * (val[xc][yc][Im_Nr] - ref_ar[i]))
                           / (float) (val[xc][yc][Im_Nr] + ref_ar[i]);
      }
      sprintf(strp, "Deviation from ref base: %.2f%%", diff_prcnt);
      txt_color("red");
      subW_TXT(loc, 5, strp);
      txt_color("blue");
   }

   x = xorigin[xc][yc] + Im_Nr*gx/(npoints-1);
   y = iHIGH - yorigin[xc][yc] - plot[xc][yc][Im_Nr];
   v_point_x = x;
   Vpointer(x, 0);
   Cpointer(x, y);
}

/* ----------------- */
   Setup_topWindow()
/* ----------------- */
{
   unsigned long white;

   if (!(XAllocNamedColor(theDisp, CMap, "white", &any_col, &rgb_col)))
   FatalError ("XAllocNamedColor problem. AJ in Setup_topWindow()");
   white = any_col.pixel;

   top_W_x = idx;
   top_W_y = GT_DLY - 2;
   topWindow = XCreateSimpleWindow(theDisp, GWindow, 0, 0,
                     top_W_x, top_W_y, 1, white, white);
   XSelectInput(theDisp, topWindow, ExposureMask);
  
}

/* --------------- */ /* Draw topWindow with text info */
   void
   DrawTopWindow()
/* --------------- */
{ 
   int  strwide;
   char *str = "Differential Image";
   char *str_cpy1 = "Lab of Brain and Cognition at NIH 1996";

   line_color("blue");
   XFillRectangle(theDisp, topWindow, theGC, 0, 0, idx, top_W_y);

   txt_color("yellow");
   XDrawString(theDisp, topWindow, txtGC, 12,
               top_W_y - 5, T_name, strlen(T_name));

   XDrawString(theDisp, topWindow, txtGC,
                  GL_DLX + 280, top_W_y - 5, str_cpy1, strlen(str_cpy1) );

   txt_color("blue");
   if ( diff_im ) {
      strwide = XTextWidth(mfinfo,str,strlen(str));
      txt_color("red");
      XDrawString(theDisp, topWindow, txtGC,
                  GL_DLX + (GX_MAX-strwide)/2, top_W_y - 5, str, strlen(str));
      txt_color("blue");
   }
}

/* ------------ */
   Setup_keys()
/* ------------ */
{
   int i;

   key = &xtkeys[0];

   for (i=0; i < N_KEYS; i++) {
       key[i].x = PADDINGW;
       key[i].y = KEY_1_Y + (i+1)*(keyhigh + PADDINGH);
       key[i].width = keywide[i];
       key[i].height = keyhigh;

       if (i > 2) {
       key[i].fore=FKeyFore;
       key[i].back=FKeyBack; }
       else {
       key[i].fore=FKeyFore;
       key[i].back=FKeyFore;  }

       key[i].wid=XCreateSimpleWindow(theDisp, GWindow, key[i].x, key[i].y,
                                      key[i].width, key[i].height, 1,
                                      key[i].fore, key[i].back);

       New_Cursor(theDisp,  key[i].wid, XC_left_ptr, "yellow", "red");
       XSelectInput(theDisp, key[i].wid, ExposureMask | ButtonPressMask |
                  ButtonReleaseMask | EnterWindowMask | LeaveWindowMask);
   }
}

/* --------------- */
   DrawKey(keynum)
   int keynum;
/* --------------- */
{
   char        *str;
   int         strwide;
   struct _key *kp;
   GC          AJkeygc; /* AJ new for proper color of the text */

   kp = &key[keynum];
   str = kp->st;
   strwide = XTextWidth(kfontinfo,str,strlen(str));

   AJkeygc=Fkeygc;

   XDrawString(theDisp ,kp->wid, AJkeygc, (kp->width-strwide)/2,
               1 + kfontinfo->ascent, str, strlen(str));
}

/* ----------------- */
   InvertKey(keynum)
   int keynum;
/* ----------------- */
{
   struct _key *kp;
   GC           AJkeyigc;

   AJkeyigc = Fkeyigc;
   kp = &key[keynum];

   XFillRectangle(theDisp, kp->wid, AJkeyigc, 0, 0, kp->width, kp->height);
}

/* ---------------- */
   LetGoKey(keynum)
   int keynum;
/* ---------------- */
{
   InvertKey(keynum);
   (*(key[keynum].fun))(keynum);
}

/* --------- */
   Ims_rot()
/* --------- */
{
   register int i, j, k, l, m, mn, n, s;

   old_im = -1;            /* to avoid faster reload of not rotated image */

   k = im_size;
   l = k - 1;
   s = ar_size;

   mn = Im_Nr*s;
   if ( rot_direct == 1 ) {                         /* clockwise rotation */
      mn = Im_Nr*s;
      for (i=0; i < k; i++) {          /* rotate and redisplay actual one */ 
         m = i * k;
         for (j=0; j < k; j++)
            a_rot[m+j] = ar[mn+(l-j)*k+i];
      }
      for (i=0; i < s; i++) ar[mn+i] = a_rot[i];
      if ( diff_im ) {                      /* rotate reference array too */
         for (i=0; i < k; i++) {
            m = i * k;
            for (j=0; j < k; j++)
               a_rot[m+j] = ref_ar[(l-j)*k+i];
         }
      }
      for (i=0; i < s; i++) ref_ar[i] = a_rot[i];
      if ( !avr_grp ) Put_image(Im_Nr);         /* put single image here */
                                                     /* rotate other data */
      for (n=0; n < N_im; n++) {
         mn = n*s;
         if ( n != Im_Nr ) {
            for (i=0; i < k; i++) {
               m = i * k;
               for (j=0; j < k; j++)
                  a_rot[m+j] = ar[mn+(l-j)*k+i];
            }
            for (i=0; i < s; i++) ar[mn+i] = a_rot[i];
         }
      }
      rot_nr += rot_direct;
      if ( avr_grp ) Put_image(Im_Nr);  /* put average image after all rot */
   }
   else if ( rot_direct == -1 ) {            /* counterclockwise rotation */
      mn = Im_Nr*s;
      for (i=0; i < k; i++) {          /* rotate and redisplay actual one */
         m = i * k;
         for (j=0; j < k; j++)
            a_rot[m+j] = ar[mn+(j+1)*k-i-1];
      }
      for (i=0; i < s; i++) ar[mn+i] = a_rot[i];
      if ( diff_im ) {                      /* rotate reference array too */
         for (i=0; i < k; i++) {
            m = i * k;
            for (j=0; j < k; j++)
               a_rot[m+j] = ref_ar[(j+1)*k-i-1];
         }
      }
      for (i=0; i < s; i++) ref_ar[i] = a_rot[i];
      if ( !avr_grp ) Put_image(Im_Nr);         /* put single image here */
                                                     /* rotate other data */
      for (n=0; n < N_im; n++) {
         mn = n*s;
         if ( n != Im_Nr ) {
            for (i=0; i < k; i++) {
               m = i * k;
               for (j=0; j < k; j++)
                  a_rot[m+j] = ar[mn+(j+1)*k-i-1];
            }
            for (i=0; i < s; i++) ar[mn+i] = a_rot[i];
         }
      }
      rot_nr += rot_direct;
      if ( avr_grp ) Put_image(Im_Nr);  /* put average image after all rot */
   }
   redraw_graph();
   DrawSubWindow();
}

/* ---------- */
   Im_diff()
/* ---------- */
{
   XUnmapWindow(theDisp,  key[kDIF].wid);
   XMapWindow  (theDisp,  key[kIR1].wid);
   XMapWindow  (theDisp,  key[kIR2].wid);
}

/* ---------- */
   Im_Aver()
/* ---------- */
{
   avr_grp = fim_avr = 0;
   XMapWindow(theDisp,  key[kAV1].wid);
   XMapWindow(theDisp,  key[kAV2].wid);
}

/* --------- */
   Im_norm()
/* --------- */
{
   diff_im = fim_dif = 0;
   avr_grp = fim_avr = 0;
   Av_length = 1;
   XUnmapWindow(theDisp,  key[kNRM].wid);
   XMapWindow  (theDisp,  key[kDIF].wid);
   old_im = -1;            /* to avoid faster reload of old data */
   Put_image(Im_Nr);
   redraw_graph();
   DrawSubWindow();
}

/* ---------- */
   Av_im1()
/* ---------- */  /* set first image for average one */
{
   Av_1 = Im_Nr;
   av1_done = 1;
}

/* ---------- */
   Av_im2()
/* ---------- */  /* set second image for average one */
{
   register int  i, j;

   if ( av1_done ) {
      Av_2 = Im_Nr;
      if ( Av_2 < Av_1 ) {
         i = Av_2; Av_2 = Av_1; Av_1 = i;
      }
      Im_Nr = Av_1;
   }
   else     Av_1 = Av_2 = Im_Nr;

   for (i=0; i< ar_size; i++) av_ar[i] = 0;
   for (i=0; i < ar_size; i++) {
      for (j=Av_1; j <= Av_2; j++) 
         av_ar[i] += ar[j*ar_size+i];
   }
   Av_length =  Av_2 - Av_1 + 1;
   for (i=0; i < ar_size; i++) av_ar[i] = av_ar[i] / Av_length;

   XUnmapWindow(theDisp,  key[kAV1].wid);
   XUnmapWindow(theDisp,  key[kAV2].wid);
   XMapWindow  (theDisp,  key[kNRM].wid);
   avr_grp = fim_avr = 1;
   av1_done = 0;
   old_im = -1;            /* to avoid faster reload of old data */
   Put_image(Im_Nr);
   redraw_graph();
   DrawSubWindow();
}

/* ---------- */
   Ref_im1()
/* ---------- */  /* set first image for average reference one */
{
   Im_1 = Im_Nr;
   im1_done = 1;
}

/* ---------- */
   Ref_im2()
/* ---------- */  /* set second image for average and make refer im for diff */
{
   register int  i, j, m;

   if ( im1_done ) {
      Im_2 = Im_Nr;
      if ( Im_2 < Im_1 ) {
         i = Im_2; Im_2 = Im_1; Im_1 = i;
      }
   }
   else     Im_1 = Im_2 = Im_Nr;

   for (i=0; i< ar_size; i++) ref_ar[i] = 0;
   for (i=0; i < ar_size; i++) {
      for (j=Im_1; j <= Im_2; j++) 
         ref_ar[i] += ar[j*ar_size+i];
   }
   m = Im_2 - Im_1 + 1;
   for (i=0; i < ar_size; i++) ref_ar[i] = ref_ar[i] / m;

   XUnmapWindow(theDisp,  key[kIR1].wid);
   XUnmapWindow(theDisp,  key[kIR2].wid);
   XMapWindow  (theDisp,  key[kNRM].wid);

   diff_im = fim_dif = 1;
   im1_done = 0;
   old_im = -1;            /* to avoid faster reload of old data */
   Put_image(Im_Nr);
   redraw_graph();
   DrawSubWindow();
   DrawTopWindow();
}

/* --------- */
   Im_help()
/* --------- */
{
   The_Help(0);
}

/* ------------ */
   draw_frame()
/* ------------ */
{
   register int i;
                               /* draw frame */
   line_color("blue");
   for (i=0;i<=mat;i++) {	   
     plotx(mdx1       ,mdy1+i*gy, 0);
     plotx(mdx1+mat*gx,mdy1+i*gy, 1);
   }
   for (i=0;i<=mat;i++) {	   
     plotx(mdx1+i*gx,mdy1, 0);
     plotx(mdx1+i*gx,mdy1+mat*gy, 1);
   }
}

/* --------------------- */	/* Reload pixmap pxWind to GWindow */
   void graphic_store()
/* --------------------- */
{
   XSetWindowBackgroundPixmap(theDisp, GWindow, pxWind);
   XClearWindow(theDisp, GWindow);
   XFlush(theDisp);
}

/* ------------------- */	/* It plots line to point (x,y) for mod = 1 */
   void plotx(x,y,mod)		/* or moves to this point for mod = 0.      */
   int x, y, mod;               /* All into the pxWind.                     */
/* ------------------- */
{
   int	iy = iHIGH - y;

   if(mod == 0) { x00 = x; y00 = iy; }
   if(mod == 1) {
     XDrawLine(theDisp, pxWind, theGC, x00, y00, x, iy);
     x00 = x;	y00 = iy;
   }
}

/* --------------------- */	/* Plot text in pxWind at x,y position */
   void plx_txt(x,y,str)        /*  relative to lower left corner (!). */
   int  x, y;
   char *str;
/* --------------------- */
{
   int	iy = iHIGH - y, n = strlen(str);;
   XDrawString(theDisp, pxWind, txtGC, x, iy, str, n);
}

/* -------------------------- */ /* Plot text in any window w at x, y */
   void plx_TXT(w, x, y, str)    /* relative to lower left corner (!)   */
   Window w;
   int  x, y;
   char *str;
/* -------------------------- */
{
   Window r;
   int x0, y0;
   u_int  width, height, bw, dp;

   if (!XGetGeometry(theDisp, w, &r, &x0, &y0, &width, &height, &bw, &dp)) {
      printf("\n Problem in plx_TXT() with XGetGeometry\n");
      exit(10);
   }
   else
      XDrawString(theDisp, w, txtGC, x, height - y, str, strlen(str));
}

/* -------------------------- */ /* Plot text in subWindow  at x, y */
   void subW_TXT(x, y, str)      /* relative to lower left corner (!) */
   int  x, y;
   char *str;
/* -------------------------- */
{
   XDrawString(theDisp, subWindow, txtGC, x, sub_W_y - y, str, strlen(str));
}

/* ----------------------- */	/* erase to background color */
   erase_graph()        /*   */
/* ----------------------- */
{
   line_color("white");
   XFillRectangle(theDisp, pxWind, theGC, 0, 0, iWIDE, iHIGH);
}

/* ----------------------- */	/* redraw entire graph */
   redraw_graph()
/* ----------------------- */
{
   erase_graph();
   draw_marker();
   draw_frame();
   plot_line();
                                      /* draw min & max values in GWindow */
   sprintf(strp, "%05d", pmax[xc][yc]); 
   plx_txt(xspace + 5, GB_DLY + GY_MAX - mytxt, strp);
   sprintf(strp, "%05d", pmin[xc][yc]);
   plx_txt(xspace + 5, GB_DLY + 5, strp);

   graphic_store();
}

/* -------------------- */     /* Change color for plotting */
   void line_color(col)        /* col - named color         */
   char *col;
/* -------------------- */
{
   XColor  any_col, rgb_col;

   if (!(XAllocNamedColor(theDisp, CMap, col, &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. ");
   XSetForeground(theDisp, theGC, any_col.pixel);
}

/* -------------------- */     /* Change color for plotting */
   void txt_color(col)         /* col - named color         */
   char *col;
/* -------------------- */
{
   XColor  any_col, rgb_col;

   if (!(XAllocNamedColor(theDisp, CMap, col, &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem for text. ");
   XSetForeground(theDisp, txtGC, any_col.pixel);
}

/* ----------------------- */
   FatalError (identifier)
   char *identifier;
/* ----------------------- */
{
   fprintf(stderr, "%s: %s\n",ProgramName, identifier);
   exit(-1);
}

/* ----------------------------------- */
   CreateGraphWindow(argv, argc)
   int  argc;
   char *argv[];
/* ----------------------------------- */
{
   XClassHint		class;
   XSetWindowAttributes attr;
   unsigned int		attrmask;
   XSizeHints		hints;
   int 			x = 0, y = 0;


   class.res_name  = "Graph";
   class.res_class = "Graph";

   hints.width = iWIDE;         hints.height = iHIGH;
   hints.max_width = iWIDE;     hints.max_height = iHIGH;
   hints.flags = PMaxSize;

   hints.min_width = iWIDE;     hints.min_height = iHIGH;
   hints.flags |= PMinSize;

   attr.background_pixel = bcol;
   attr.border_pixel     = fcol;
   attrmask = CWBackPixel | CWBorderPixel;

   GWindow = XCreateWindow(theDisp, rootW, x, y, iWIDE, iHIGH, 2,
	CopyFromParent, CopyFromParent, CopyFromParent, attrmask, &attr);

   if (!GWindow)
      FatalError("unable to open window (X11 windows not running).");

   XSetClassHint(theDisp, GWindow, &class);
   XSetStandardProperties(theDisp, GWindow, "Time_Course", "Time_Course", None,
			  argv, argc, &hints);
}

/* -------------- */     /* discard events x (of ev) to stop faster */
   discard(x, ev)
   int x;
   XEvent *ev;
/* -------------- */
{
   while ( XCheckWindowEvent(theDisp, theWindow, x, ev) ) ;
   while ( XCheckWindowEvent(theDisp,   GWindow, x, ev) ) ;
}

/* ------------------------ */
   discard_Key(keyW, x, ev)
   int    x;
   XEvent *ev;
   Window keyW;
/* ------------------------ */
{
   while ( XCheckWindowEvent(theDisp, keyW, x, ev) ) ;
}

/* -------------------- */
   Track_Cursor(mx, my)
   int mx, my;
/* -------------------- */
{
   Window       rW, cW;
   u_int        key;
   int          x, y, rx, ry;

   while (XQueryPointer(theDisp, theWindow, &rW, &cW,
                                &rx, &ry, &x, &y, &key)) {
      if ( !(key & Button1Mask) ) break;    /* button released */

      if ( mx != x || my != y ) {   /* this marker was moved */
         xpoint = Mltx[x]/x_mag;
         ypoint = Mlty[y]/x_mag;
         DrawSubWindow();
      }
      mx = x;
      my = y;
   }
}
/* -------------- */
   Vpointer(x, y)
   int x, y;
/* -------------- */
{
   XPoint a[3];
                                 /* fill pattern is shifted by 1 pixel => -1 */
   a[0].x = x-1;    a[0].y = y-1;
   a[1].x = x-8;    a[1].y = y+14;
   a[2].x = x+6;    a[2].y = y+14;

   line_color("red");
   XFillPolygon(theDisp, subWindow, theGC, a, 3, Convex, CoordModeOrigin);
}

/* -------------- */
   Cpointer(x, y)
   int x, y;
/* -------------- */
{
   int  i;
   XPoint a[12];

   for (i=0; i < 12; i++) {
      a[i].x = sm_cir[i].x + x;
      a[i].y = sm_cir[i].y + y;
   }
     
   line_color("red");
   XDrawPoints(theDisp, GWindow, theGC, a, 12, CoordModeOrigin);
}

/* ---------------- */ /* tracks continuousely v pointer and sets image nr */
   Track_Vpointer()
/* ---------------- */
{
   int          Im_Old, im, c_im;
   Window       rW, cW;
   u_int        key;
   int          x, y, rx, ry;
   int          aaa = avr_grp, ddd = diff_im, redr;

   im = Im_Old = Im_Nr;

   New_Cursor(theDisp, subWindow, XC_left_ptr, "red", "white");

   while (XQueryPointer(theDisp, subWindow, &rW, &cW,
                                &rx, &ry, &x, &y, &key)) {
      if ( !(key & Button1Mask) ) break;    /* button released */

      if ( v_point_x != x ) {   /* the marker was moved */

         c_im = (min(gx, max(0, x - xorigin[xc][yc]))*npoints) / gx;
         if ( c_im < 0 )       c_im = 0;
         if ( c_im > npoints - Av_length ) c_im = npoints - Av_length;
         if ( im != c_im ) {
            Im_Nr = c_im;
            XClearWindow(theDisp, GWindow);
            DrawSubWindow();
         }
         im = c_im;
      }
   }
   v_point_x = xorigin[xc][yc] + Im_Nr*gx/(npoints-1);

   if ( Im_Nr != Im_Old ) {
      diff_im = fim_dif;
      avr_grp = fim_avr;
      Put_image(Im_Nr);
      redr =  avr_grp * 4 + (aaa - avr_grp) *2 + ddd - diff_im;
      if ( redr ) {
         old_im = -1;    /* to avoid faster reload of old data */
         redraw_graph();
         DrawTopWindow();
         DrawSubWindow();
      }
   }

   New_Cursor(theDisp, subWindow, XC_left_ptr, "blue", "yellow");
}

/* ------------------ */
   int is_file(fname)
   char  *fname;
/* ------------------ */
{
   FILE          *fp;

   if ( (fp = fopen(fname, "r")) != NULL ) {    /* return = 1 if file exist */
      fclose(fp);
      return(1);
   }
   else
      return(0);
}

/* ----------------------------- */   /* write plot to file  */
   void print_plot()
/* ----------------------------- */
{
   int  i, x, y;

   if ( txtW_ON ) {
      XBell(theDisp, 100);
      return;
   }
   txtW_ON = 1;

   x = 50 + GL_DLX;
   y = 50 + GT_DLY;

   strp[0] = NULL;
   take_file_name(theDisp, GWindow , CMap, txtGC, mfinfo, x, y, strp, 41,
                  "Enter plot name:");

   if ( strp[0] != NULL ) {
      sprintf(plotbuf,"%d\n", val[xc][yc][0]);
      for (i=1; i < npoints; i++) {
         sprintf(fnum,"%d\n", val[xc][yc][i]);
	     strcat(plotbuf, fnum);
      }
   
      isize = strlen(plotbuf);
      i = WRite_iqm(strp, &isize, plotbuf);
      if ( i != 0 ) XBell(theDisp, 100);
   }
   txtW_ON = 0;
}

/* ----------------- */
   int save_act_im()
/* ----------------- */
{
   int  i, x, y;
   XWindowAttributes  wat;
   Window             ww;

   if ( txtW_ON ) {
      XBell(theDisp, 100);
      return(2);
   }
   txtW_ON = 1;

/*   OK in relation to the root window
   XGetWindowAttributes(theDisp, theWindow, &wat);
   XTranslateCoordinates(theDisp, theWindow, wat.root, -wat.border_width,
                        -wat.border_width, &x, &y, &ww);
   x += eWIDE/2;
   y += eHIGH/2;

   strp[0] = NULL;
   take_file_name(theDisp, rootW, CMap, txtGC, mfinfo, x, y, strp, 40,
                  "Enter image name:");
*/

   x = 50 + GL_DLX;
   y = 50 + GT_DLY;

   strp[0] = NULL;
   take_file_name(theDisp, GWindow , CMap, txtGC, mfinfo, x, y, strp, 41,
                  "Enter image name:");
   if ( strp[0] != NULL ) {
      i = WRite_iqm(strp, &fsize, tmp_ar);
      if ( i != 0 ) XBell(theDisp, 100);
   }

   txtW_ON = 0;
   return(0);
}

/* popup window which take a file name of max length str_l */
/* --------------------------------------- return name of file in name */
   int take_file_name(theDisp, topW, CMap, txtGC, finf, x, y, name, str_l, text)
   Display     *theDisp;
   Window      topW;                    /* top window of this popup one */
   Colormap    CMap;
   GC          txtGC;
   XFontStruct *finf;
   int         x, y;                    /* relative position */
   char *name, *text;                   /* name - file name, max str_l chars */
   int  str_l;                          /* text - window's header */
/* --------------------------------------- */
{
   int    w_l1, w_h1, w_l2, w_h2, h_txt, expose = 0, name_OK = 0, length;
   int    x_txt, y_txt, x_l2, y_l2, error = 0;
   int    eee[2];
   char   *errr = "File exist:", *DT = "_";

   XEvent ev;
   XColor             any_col, rgb_col;
   unsigned long      border, back1, back2;
   Window             txtWindow1, txtWindow2, ww;
   XSizeHints         hints;

   eee[0] = eee[1] = 0;

   if (!(XAllocNamedColor(theDisp, CMap, "red", &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. AJ in save_act_im()");
   border = any_col.pixel;
   if (!(XAllocNamedColor(theDisp, CMap, "yellow", &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. AJ in save_act_im()");
   back1= any_col.pixel;
   if (!(XAllocNamedColor(theDisp, CMap, "white", &any_col, &rgb_col)))
      FatalError ("XAllocNamedColor problem. AJ in save_act_im()");
   back2= any_col.pixel;

   hints.flags = USPosition;
   hints.x = x;
   hints.y = y;
   
   h_txt =  finf->max_bounds.ascent + finf->max_bounds.descent;
   w_l1 = (str_l - 1) * finf->max_bounds.width + 30;
   w_h1 = 3 * h_txt + 18;
   x_l2 = 10;
   y_l2 = 2 * h_txt + 5;
   w_l2 = w_l1 - 20;
   w_h2 = h_txt + 8;
   x_txt = 5;
   y_txt = h_txt;

   txtWindow1 =  XCreateSimpleWindow(theDisp, topW, x, y, w_l1, w_h1, 1,
                                     border, back1);
   XSelectInput(theDisp, txtWindow1, ExposureMask);
   XSetStandardProperties(theDisp, txtWindow1 , "Text", "Text", None,
                          NULL, NULL, &hints);
   txtWindow2 =  XCreateSimpleWindow(theDisp, txtWindow1, x_l2, y_l2,
                            w_l2, w_h2, 1, border, back2);
   XSelectInput(theDisp, txtWindow2, ExposureMask | KeyPressMask);
   New_Cursor(theDisp, txtWindow1, XC_left_ptr, "blue", "white");
   New_Cursor(theDisp, txtWindow2, XC_xterm, "blue", "white");
   XMapWindow(theDisp, txtWindow1);
   XMapWindow(theDisp, txtWindow2);

   while ( !expose ) {                /* wait for expose */
      Window wind;
      XNextEvent(theDisp, &ev);
      switch (ev.type) {
         case Expose: {
            XExposeEvent *e = (XExposeEvent *) &ev;
            wind = e->window;
            if (wind ==  txtWindow1) eee[0] = 1;
            if (wind ==  txtWindow2) eee[1] = 1;
            expose = eee[0] * eee[1];
         }
         default:                /* ignore unexpected events */
            break;
      }
   } /* end of while( !expose ) */

   txt_color("red");
   XDrawString(theDisp, txtWindow1, txtGC, 25, h_txt + 3 ,
               text, strlen(text));
   strncat(name, DT, 1);
   txt_color("blue");
   XDrawString(theDisp, txtWindow2, txtGC, x_txt, y_txt,
               name, strlen(name));

   txt_color("blue");

   while ( !name_OK ) {                /* wait for file name */
      Window wind;
      XNextEvent(theDisp, &ev);
      switch (ev.type) {
         case Expose: {
            XExposeEvent *e = (XExposeEvent *) &ev;

            wind = e->window;
            if (wind ==  txtWindow1) {
               txt_color("red");
               XDrawString(theDisp, txtWindow1, txtGC, 30, h_txt + 3 ,
                           text, strlen(text));
               txt_color("blue");
            }
            if (wind ==  txtWindow2) {
               txt_color("blue");
               XDrawString(theDisp, txtWindow2, txtGC, x_txt, y_txt,
                           name, strlen(name));
               txt_color("blue");
            }
         }
         case KeyPress: {
            XKeyEvent *key_event = (XKeyEvent *) &ev;
            u_char      buf[128];
            KeySym    ks;
            XComposeStatus status;

            wind = key_event->window;
            if (wind ==  txtWindow2) {
               if (error) {
                  XClearWindow(theDisp, txtWindow1);
                  txt_color("red");
                  XDrawString(theDisp, txtWindow1, txtGC, 25, h_txt + 3 ,
                              text, strlen(text));
                  txt_color("blue");
                  error = 0;
               }
               buf[0] = 0;
               XLookupString(key_event, buf, 128, &ks, &status);
               if ( (ks == XK_Return)   ||
                    (ks == XK_Linefeed) ) {
                  length = strlen(name);
                  name[length-1] = NULL;
                  if ( name[0] == NULL ) name_OK = 1;
                  else if ( is_file(name) ) {
                     XClearWindow(theDisp, txtWindow1);
                     txt_color("red");
                     XDrawString(theDisp, txtWindow1, txtGC, 10, h_txt ,
                           errr, strlen(errr));
                     XDrawString(theDisp, txtWindow1, txtGC, 15, 2*h_txt ,
                           name, strlen(name));
                     txt_color("blue");
                     name[0] = NULL;
                     strncat(name, DT, 1);
                     XBell(theDisp, 100);
                     error = 1;
                  }
                  else
                     name_OK = 1;
               }
               else if ( ((ks >= XK_plus) && (ks <= XK_9)) ||
                         ((ks >= XK_A) && (ks <= XK_Z))  ||
                         ((ks >= XK_asciicircum) && (ks <= XK_z)) ||
                          (ks == XK_asciitilde) ) {
                  if ( (length = strlen(name)) + 1 > str_l ) {
                     name[length-1] = ' ';
                     XBell(theDisp, 100);
                  }
                  else if ( ((buf[0] >= 43) && (buf[0] <= 57)) ||
                            ((buf[0] >= 65) && (buf[0] <= 90)) ||
                            ((buf[0] >= 94) && (buf[0] <= 122))||
                             (buf[0] == 126) ) {
                     name[length - 1] = NULL;
                     strncat(name, (char *) buf, 1);
                     strncat(name, DT, 1);
                  }
               }
               else if ( ((ks >= XK_Shift_L) && (ks <= XK_Hyper_R)) ||
                         ((ks >= XK_F1) && (ks <= XK_F35)) ||
                         ((ks >= XK_KP_0) && (ks <= XK_KP_9)) )
                     ;     /* do nothing */
               else if ((ks == XK_BackSpace) || (ks == XK_Delete)) {
                  if ( (length = strlen(name)) > 1 ) {
                     name[length - 2] = NULL;
                     strncat(name, DT, 1);
                  }
                  else
                     XBell(theDisp, 100);
               }
               else
                  XBell(theDisp, 100);
            }
            XClearWindow(theDisp, txtWindow2);
            txt_color("blue");
            XDrawString(theDisp, txtWindow2, txtGC, x_txt, y_txt,
                        name, strlen(name));
            txt_color("blue");
         }
         default:                /* ignore unexpected events */
            break;
      }

   } /* end of while( !name_OK ) */

   txt_color("blue");
   XDestroyWindow(theDisp, txtWindow1);
}
