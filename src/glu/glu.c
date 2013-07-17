#include <GLES/gl.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

void gluPerspective(double fovy, double aspect, double near, double far) {
    printf("gluPerspective(%.1f, %.1f, %.1f, %.1f);\n", fovy, aspect, near, far);
    GLfloat m[16];
    double si, co, dz;
    double rad = fovy / 2 * M_PI / 180;
    double a, b, c, d;

    dz = far - near;
    si = sin(rad);
    if (dz == 0 || si == 0 || aspect == 0)
        return;
    co = cos(rad) / si;

    a = co / aspect;
    b = co;
    c = -(far + near) / dz;
    d = -2 * near * far / dz;

    # define M(X,Y)  m[Y * 4 + X]
    M(0,0) = a; M(0,1) = 0; M(0,2) = 0;  M(0,3) = 0;
    M(1,0) = 0; M(1,1) = b; M(1,2) = 0;  M(1,3) = 0;
    M(2,0) = 0; M(2,1) = 0; M(2,2) = c;  M(2,3) = d;
    M(3,0) = 0; M(3,1) = 0; M(3,2) = -1; M(3,3) = 0;
    # undef M

    glMultMatrixf (m);
}

GLint gluBuild2DMipmaps (GLenum target, GLint internalFormat, GLsizei width, GLsizei height, GLenum format, GLenum type, const void *data) {
  return 0;
}

void gluOrtho2D (double left, double right, double bottom, double top) {

    GLfloat m[16];

    double a, b, c, d;

    a = 2.0 / (right - left);
    b = 2.0 / (top - bottom);
    c = - (right + left) / (right - left);
    d = - (top + bottom) / (top - bottom);

    # define M(X,Y)  m[Y * 4 + X]
    M(0,0) = a; M(0,1) = 0; M(0,2) = 0;  M(0,3) = 0;
    M(1,0) = 0; M(1,1) = b; M(1,2) = 0;  M(1,3) = 0;
    M(2,0) = 0; M(2,1) = 0; M(2,2) = -1;  M(2,3) = 0;
    M(3,0) = c; M(3,1) = d; M(3,2) = 0; M(3,3) = 1;
    # undef M

    glMultMatrixf (m);
}

void gluLookAt( double eyex, double eyey, double eyez, 
                double centerx, double centery, double centerz, 
                double upx, double upy, double upz ) 
{ 
   double m[16]; 
   double x[3], y[3], z[3]; 
   double mag; 
   /* Make rotation matrix */ 
   /* Z vector */ 
   z[0] = eyex - centerx; 
   z[1] = eyey - centery; 
   z[2] = eyez - centerz; 
   mag = sqrt( z[0]*z[0] + z[1]*z[1] + z[2]*z[2] ); 
   if (mag) {  /* mpichler, 19950515 */ 
      z[0] /= mag; 
      z[1] /= mag; 
      z[2] /= mag; 
   } 
   /* Y vector */ 
   y[0] = upx; 
   y[1] = upy; 
   y[2] = upz; 
   /* X vector = Y cross Z */ 
   x[0] =  y[1]*z[2] - y[2]*z[1]; 
   x[1] = -y[0]*z[2] + y[2]*z[0]; 
   x[2] =  y[0]*z[1] - y[1]*z[0]; 
   /* Recompute Y = Z cross X */ 
   y[0] =  z[1]*x[2] - z[2]*x[1]; 
   y[1] = -z[0]*x[2] + z[2]*x[0]; 
   y[2] =  z[0]*x[1] - z[1]*x[0]; 
   /* mpichler, 19950515 */ 
   /* cross product gives area of parallelogram, which is < 1.0 for 
    * non-perpendicular unit-length vectors; so normalize x, y here 
    */ 
   mag = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ); 
   if (mag) { 
      x[0] /= mag; 
      x[1] /= mag; 
      x[2] /= mag; 
   } 
   mag = sqrt( y[0]*y[0] + y[1]*y[1] + y[2]*y[2] ); 
   if (mag) { 
      y[0] /= mag; 
      y[1] /= mag; 
      y[2] /= mag; 
   } 
#define M(row,col)  m[col*4+row] 
   M(0,0) = x[0];  M(0,1) = x[1];  M(0,2) = x[2];  M(0,3) = 0.0; 
   M(1,0) = y[0];  M(1,1) = y[1];  M(1,2) = y[2];  M(1,3) = 0.0; 
   M(2,0) = z[0];  M(2,1) = z[1];  M(2,2) = z[2];  M(2,3) = 0.0; 
   M(3,0) = 0.0;   M(3,1) = 0.0;   M(3,2) = 0.0;   M(3,3) = 1.0; 
#undef M 
   glMultMatrixd( m ); 
   /* Translate Eye to Origin */ 
   glTranslated( -eyex, -eyey, -eyez );  
}
