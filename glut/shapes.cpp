#include "screencasts.h"



/*
*  cube
*  ------
*  Draw a cube
*     at (x,y,z)
*     dimensions (dx,dy,dz)
*     rotated th about the y axis
*/
void cube(double x, double y, double z,double dx, double dy, double dz,double th,int i,int j,int k)
{
	if (dronpa[i+1][j+1][k+1].ter == true)  glColor4f(0, 1, 0, 1);	//zöld szín ha térrel terhelt molekulát adunk meg
	else glColor4f(1, 0, 0, 1);		//alap piros szín

	//nagyság tényezõ
	double s = 0.5;
	/*  Cube vertices */
	GLfloat A[3] = { 0.5 + s / 2, 0.5+s/2 , 0.5 + s / 2 };
	GLfloat B[3] = { 0.5 - s / 2, 0.5+s/2 , 0.5 + s / 2 };
	GLfloat C[3] = { 0.5 - s / 2,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat D[3] = { 0.5 + s / 2,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat E[3] = { 0.5 + s / 2, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat F[3] = { 0.5 - s / 2, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat G[3] = { 0.5 - s / 2,0.5 - s / 2,0.5 - s / 2 };
	GLfloat H[3] = { 0.5 + s / 2,0.5 - s / 2,0.5 - s / 2 };

	GLfloat AX[3] = { 0.5 + s / 2-s/3, 0.5 + s / 2 , 0.5 + s / 2 };
	GLfloat BX[3] = { 0.5 - s / 2+s/3, 0.5 + s / 2 , 0.5 + s / 2 };
	GLfloat CX[3] = { 0.5 - s / 2+s/3,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat DX[3] = { 0.5 + s / 2-s/3,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat EX[3] = { 0.5 + s / 2-s/3, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat FX[3] = { 0.5 - s / 2+s/3, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat GX[3] = { 0.5 - s / 2+s/3,0.5 - s / 2,0.5 - s / 2 };
	GLfloat HX[3] = { 0.5 + s / 2-s/3,0.5 - s / 2,0.5 - s / 2 };

	GLfloat AZ[3] = { 0.5 + s / 2, 0.5 + s / 2 , 0.5 + s / 2-s/3 };
	GLfloat BZ[3] = { 0.5 - s / 2, 0.5 + s / 2 , 0.5 + s / 2-s/3 };
	GLfloat CZ[3] = { 0.5 - s / 2,0.5 - s / 2, 0.5 + s / 2-s/3 };
	GLfloat DZ[3] = { 0.5 + s / 2,0.5 - s / 2, 0.5 + s / 2-s/3 };
	GLfloat EZ[3] = { 0.5 + s / 2, 0.5 + s / 2,0.5 - s / 2+s/3 };
	GLfloat FZ[3] = { 0.5 - s / 2, 0.5 + s / 2,0.5 - s / 2+s/3 };
	GLfloat GZ[3] = { 0.5 - s / 2,0.5 - s / 2,0.5 - s / 2+s/3 };
	GLfloat HZ[3] = { 0.5 + s / 2,0.5 - s / 2,0.5 - s / 2+s/3 };

	
	glPushMatrix();
	/*  Transform */
	glTranslated(x, y, z);
	glRotated(th, 0, 1, 0);
	glScaled(dx, dy, dz);

	glEnable(GL_TEXTURE_2D);
	/* using the current texture */
	

	/* Cube */
	/* front => ABCD */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(C);
	glTexCoord2f(1, 0); glVertex3fv(CX);
	glTexCoord2f(1, 1); glVertex3fv(BX);
	glTexCoord2f(0, 1); glVertex3fv(B);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(CX);
	glTexCoord2f(1, 0); glVertex3fv(DX);
	glTexCoord2f(1, 1); glVertex3fv(AX);
	glTexCoord2f(0, 1); glVertex3fv(BX);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(DX);
	glTexCoord2f(1, 0); glVertex3fv(D);
	glTexCoord2f(1, 1); glVertex3fv(A);
	glTexCoord2f(0, 1); glVertex3fv(AX);
	glEnd();

	/* back => FEHG */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(F);
	glTexCoord2f(1, 0); glVertex3fv(FX);
	glTexCoord2f(1, 1); glVertex3fv(GX);
	glTexCoord2f(0, 1); glVertex3fv(G);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(FX);
	glTexCoord2f(1, 0); glVertex3fv(EX);
	glTexCoord2f(1, 1); glVertex3fv(HX);
	glTexCoord2f(0, 1); glVertex3fv(GX);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(EX);
	glTexCoord2f(1, 0); glVertex3fv(E);
	glTexCoord2f(1, 1); glVertex3fv(H);
	glTexCoord2f(0, 1); glVertex3fv(HX);
	glEnd();

	/* right => EADH */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(D);
	glTexCoord2f(1, 0); glVertex3fv(DZ);
	glTexCoord2f(1, 1); glVertex3fv(AZ);
	glTexCoord2f(0, 1); glVertex3fv(A);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(DZ);
	glTexCoord2f(1, 0); glVertex3fv(HZ);
	glTexCoord2f(1, 1); glVertex3fv(EZ);
	glTexCoord2f(0, 1); glVertex3fv(AZ);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(HZ);
	glTexCoord2f(1, 0); glVertex3fv(H);
	glTexCoord2f(1, 1); glVertex3fv(E);
	glTexCoord2f(0, 1); glVertex3fv(EZ);
	glEnd();

	/* left => BFGC */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(G);
	glTexCoord2f(1, 0); glVertex3fv(GZ);
	glTexCoord2f(1, 1); glVertex3fv(FZ);
	glTexCoord2f(0, 1); glVertex3fv(F);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(GZ);
	glTexCoord2f(1, 0); glVertex3fv(CZ);
	glTexCoord2f(1, 1); glVertex3fv(BZ);
	glTexCoord2f(0, 1); glVertex3fv(FZ);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(CZ);
	glTexCoord2f(1, 0); glVertex3fv(C);
	glTexCoord2f(1, 1); glVertex3fv(B);
	glTexCoord2f(0, 1); glVertex3fv(BZ);
	glEnd();

	/* top => EFBA */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(B);
	glTexCoord2f(1, 0); glVertex3fv(BX);
	glTexCoord2f(1, 1); glVertex3fv(FX);
	glTexCoord2f(0, 1); glVertex3fv(F);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(BX);
	glTexCoord2f(1, 0); glVertex3fv(AX);
	glTexCoord2f(1, 1); glVertex3fv(EX);
	glTexCoord2f(0, 1); glVertex3fv(FX);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(AX);
	glTexCoord2f(1, 0); glVertex3fv(A);
	glTexCoord2f(1, 1); glVertex3fv(E);
	glTexCoord2f(0, 1); glVertex3fv(EX);
	glEnd();

	/* bottom => DCGH */
	currentTexture = textures[i];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(G);
	glTexCoord2f(1, 0); glVertex3fv(GX);
	glTexCoord2f(1, 1); glVertex3fv(CX);
	glTexCoord2f(0, 1); glVertex3fv(C);
	glEnd();

	currentTexture = textures[j];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(GX);
	glTexCoord2f(1, 0); glVertex3fv(HX);
	glTexCoord2f(1, 1); glVertex3fv(DX);
	glTexCoord2f(0, 1); glVertex3fv(CX);
	glEnd();

	currentTexture = textures[k];
	glBindTexture(GL_TEXTURE_2D, currentTexture);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glTexCoord2f(0, 0); glVertex3fv(HX);
	glTexCoord2f(1, 0); glVertex3fv(H);
	glTexCoord2f(1, 1); glVertex3fv(D);
	glTexCoord2f(0, 1); glVertex3fv(DX);
	glEnd();

	glPopMatrix();
	glDisable(GL_TEXTURE_2D);
	glColor4f(1, 1, 1, 1);
}


//ugyanaz mint a cube függvény, csak kisebb fekete kockákat rajzol
void negyzetracs(double x, double y, double z,double dx, double dy, double dz,double th,int i, int j,int k)
{
	
	double s = 0.2;
	/*  Cube vertices */
	GLfloat A[3] = { 0.5 + s / 2  , 0.5 + s / 2   , 0.5 + s / 2   };
	GLfloat B[3] = { 0.5 - s / 2  , 0.5 + s / 2  , 0.5 + s / 2   };
	GLfloat C[3] = { 0.5 - s / 2  ,0.5 - s / 2  , 0.5 + s / 2  };
	GLfloat D[3] = { 0.5 + s / 2 ,0.5 - s / 2  , 0.5 + s / 2  };
	GLfloat E[3] = { 0.5 + s / 2  , 0.5 + s / 2  ,0.5 - s / 2   };
	GLfloat F[3] = { 0.5 - s / 2  , 0.5 + s / 2  ,0.5 - s / 2   };
	GLfloat G[3] = { 0.5 - s / 2 ,0.5 - s / 2  ,0.5 - s / 2   };
	GLfloat H[3] = { 0.5 + s / 2  ,0.5 - s / 2  ,0.5 - s / 2   };

	GLfloat AX[3] = { 0.5 + s / 2 - s / 3, 0.5 + s / 2 , 0.5 + s / 2 };
	GLfloat BX[3] = { 0.5 - s / 2 + s / 3, 0.5 + s / 2 , 0.5 + s / 2 };
	GLfloat CX[3] = { 0.5 - s / 2 + s / 3,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat DX[3] = { 0.5 + s / 2 - s / 3,0.5 - s / 2, 0.5 + s / 2 };
	GLfloat EX[3] = { 0.5 + s / 2 - s / 3, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat FX[3] = { 0.5 - s / 2 + s / 3, 0.5 + s / 2,0.5 - s / 2 };
	GLfloat GX[3] = { 0.5 - s / 2 + s / 3,0.5 - s / 2,0.5 - s / 2 };
	GLfloat HX[3] = { 0.5 + s / 2 - s / 3,0.5 - s / 2,0.5 - s / 2 };

	GLfloat AZ[3] = { 0.5 + s / 2, 0.5 + s / 2 , 0.5 + s / 2 - s / 3 };
	GLfloat BZ[3] = { 0.5 - s / 2, 0.5 + s / 2 , 0.5 + s / 2 - s / 3 };
	GLfloat CZ[3] = { 0.5 - s / 2,0.5 - s / 2, 0.5 + s / 2 - s / 3 };
	GLfloat DZ[3] = { 0.5 + s / 2,0.5 - s / 2, 0.5 + s / 2 - s / 3 };
	GLfloat EZ[3] = { 0.5 + s / 2, 0.5 + s / 2,0.5 - s / 2 + s / 3 };
	GLfloat FZ[3] = { 0.5 - s / 2, 0.5 + s / 2,0.5 - s / 2 + s / 3 };
	GLfloat GZ[3] = { 0.5 - s / 2,0.5 - s / 2,0.5 - s / 2 + s / 3 };
	GLfloat HZ[3] = { 0.5 + s / 2,0.5 - s / 2,0.5 - s / 2 + s / 3 };
	
	glRotated(th, 0, 1, 0);
	glScaled(dx, dy, dz);
	glEnable(GL_TEXTURE_2D);
	glNormal3f(0, 0, 1);

	
	
	//36x36x36os rács
	for (int tex = 0; tex < 36; tex++)
	{
		currentTexture = textures[tex];
		glBindTexture(GL_TEXTURE_2D, currentTexture);
		glBegin(GL_QUADS);
		int hely = 0;
		for (int kk = -18; kk < 18; kk++)
		{
			for (int jj = -18; jj < 18; jj++)
			{
				for (int ii = -18; ii < 18; ii++)
				{
					//glPushMatrix();

					/*  Transform */
					//glTranslated(ii,jj, kk);

					if (tex == h1[hely])
					{
						/*  Cube vertices */
						GLfloat A[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat B[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat C[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat D[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat E[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat F[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat G[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat H[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat BX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat CX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat DX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat EX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat FX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat GX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat HX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat BZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat CZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat DZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat EZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat FZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat GZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat HZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };

						/* Cube */
						/* front => ABCD */
						/* back => FEHG */
						/* right => EADH */
						/* left => BFGC */
						/* top => EFBA */
						/* bottom => DCGH */


						glTexCoord2f(0, 0); glVertex3fv(C);
						glTexCoord2f(1, 0); glVertex3fv(CX);
						glTexCoord2f(1, 1); glVertex3fv(BX);
						glTexCoord2f(0, 1); glVertex3fv(B);

						glTexCoord2f(0, 0); glVertex3fv(F);
						glTexCoord2f(1, 0); glVertex3fv(FX);
						glTexCoord2f(1, 1); glVertex3fv(GX);
						glTexCoord2f(0, 1); glVertex3fv(G);

						glTexCoord2f(0, 0); glVertex3fv(D);
						glTexCoord2f(1, 0); glVertex3fv(DZ);
						glTexCoord2f(1, 1); glVertex3fv(AZ);
						glTexCoord2f(0, 1); glVertex3fv(A);

						glTexCoord2f(0, 0); glVertex3fv(G);
						glTexCoord2f(1, 0); glVertex3fv(GZ);
						glTexCoord2f(1, 1); glVertex3fv(FZ);
						glTexCoord2f(0, 1); glVertex3fv(F);

						glTexCoord2f(0, 0); glVertex3fv(B);
						glTexCoord2f(1, 0); glVertex3fv(BX);
						glTexCoord2f(1, 1); glVertex3fv(FX);
						glTexCoord2f(0, 1); glVertex3fv(F);

						glTexCoord2f(0, 0); glVertex3fv(G);
						glTexCoord2f(1, 0); glVertex3fv(GX);
						glTexCoord2f(1, 1); glVertex3fv(CX);
						glTexCoord2f(0, 1); glVertex3fv(C);
					}

					if (tex == h2[hely])
					{
						/*  Cube vertices */
						GLfloat A[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat B[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat C[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat D[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat E[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat F[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat G[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat H[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat BX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat CX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat DX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat EX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat FX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat GX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat HX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat BZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat CZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat DZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat EZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat FZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat GZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat HZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };

						/* Cube */
						/* front => ABCD */
						/* back => FEHG */
						/* right => EADH */
						/* left => BFGC */
						/* top => EFBA */
						/* bottom => DCGH */
						glTexCoord2f(0, 0); glVertex3fv(CX);
						glTexCoord2f(1, 0); glVertex3fv(DX);
						glTexCoord2f(1, 1); glVertex3fv(AX);
						glTexCoord2f(0, 1); glVertex3fv(BX);

						glTexCoord2f(0, 0); glVertex3fv(FX);
						glTexCoord2f(1, 0); glVertex3fv(EX);
						glTexCoord2f(1, 1); glVertex3fv(HX);
						glTexCoord2f(0, 1); glVertex3fv(GX);

						glTexCoord2f(0, 0); glVertex3fv(DZ);
						glTexCoord2f(1, 0); glVertex3fv(HZ);
						glTexCoord2f(1, 1); glVertex3fv(EZ);
						glTexCoord2f(0, 1); glVertex3fv(AZ);

						glTexCoord2f(0, 0); glVertex3fv(GZ);
						glTexCoord2f(1, 0); glVertex3fv(CZ);
						glTexCoord2f(1, 1); glVertex3fv(BZ);
						glTexCoord2f(0, 1); glVertex3fv(FZ);

						glTexCoord2f(0, 0); glVertex3fv(BX);
						glTexCoord2f(1, 0); glVertex3fv(AX);
						glTexCoord2f(1, 1); glVertex3fv(EX);
						glTexCoord2f(0, 1); glVertex3fv(FX);

						glTexCoord2f(0, 0); glVertex3fv(GX);
						glTexCoord2f(1, 0); glVertex3fv(HX);
						glTexCoord2f(1, 1); glVertex3fv(DX);
						glTexCoord2f(0, 1); glVertex3fv(CX);
					}

					if (tex == h3[hely])
					{
						/*  Cube vertices */
						GLfloat A[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat B[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat C[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat D[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat E[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat F[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat G[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat H[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 + kk };
						GLfloat BX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat CX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat DX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 + kk };
						GLfloat EX[3] = { 0.5 + s / 2 - s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat FX[3] = { 0.5 - s / 2 + s / 3 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat GX[3] = { 0.5 - s / 2 + s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };
						GLfloat HX[3] = { 0.5 + s / 2 - s / 3 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + kk };

						GLfloat AZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat BZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj , 0.5 + s / 2 - s / 3 + kk };
						GLfloat CZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat DZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj, 0.5 + s / 2 - s / 3 + kk };
						GLfloat EZ[3] = { 0.5 + s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat FZ[3] = { 0.5 - s / 2 + ii, 0.5 + s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat GZ[3] = { 0.5 - s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };
						GLfloat HZ[3] = { 0.5 + s / 2 + ii,0.5 - s / 2 + jj,0.5 - s / 2 + s / 3 + kk };

						/* Cube */
						/* front => ABCD */
						/* back => FEHG */
						/* right => EADH */
						/* left => BFGC */
						/* top => EFBA */
						/* bottom => DCGH */
						glTexCoord2f(0, 0); glVertex3fv(DX);
						glTexCoord2f(1, 0); glVertex3fv(D);
						glTexCoord2f(1, 1); glVertex3fv(A);
						glTexCoord2f(0, 1); glVertex3fv(AX);

						glTexCoord2f(0, 0); glVertex3fv(EX);
						glTexCoord2f(1, 0); glVertex3fv(E);
						glTexCoord2f(1, 1); glVertex3fv(H);
						glTexCoord2f(0, 1); glVertex3fv(HX);

						glTexCoord2f(0, 0); glVertex3fv(HZ);
						glTexCoord2f(1, 0); glVertex3fv(H);
						glTexCoord2f(1, 1); glVertex3fv(E);
						glTexCoord2f(0, 1); glVertex3fv(EZ);

						glTexCoord2f(0, 0); glVertex3fv(CZ);
						glTexCoord2f(1, 0); glVertex3fv(C);
						glTexCoord2f(1, 1); glVertex3fv(B);
						glTexCoord2f(0, 1); glVertex3fv(BZ);

						glTexCoord2f(0, 0); glVertex3fv(AX);
						glTexCoord2f(1, 0); glVertex3fv(A);
						glTexCoord2f(1, 1); glVertex3fv(E);
						glTexCoord2f(0, 1); glVertex3fv(EX);

						glTexCoord2f(0, 0); glVertex3fv(HX);
						glTexCoord2f(1, 0); glVertex3fv(H);
						glTexCoord2f(1, 1); glVertex3fv(D);
						glTexCoord2f(0, 1); glVertex3fv(DX);
					}
						

					//glPopMatrix();
					
					hely++;
				}
			}
		}
		glEnd();
	}
	
	
	glDisable(GL_TEXTURE_2D);
}