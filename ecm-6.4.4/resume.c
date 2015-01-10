/* Functions for reading a writing resume file lines.

Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2010, 2011, 2012
Paul Zimmermann, Alexander Kruppa and Cyril Bouvier.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include <gmp.h>
#include "ecm.h"
#include "ecm-ecm.h"
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#if defined (_MSC_VER) || defined (__MINGW32__)
/* needed to declare GetComputerName() for write_resumefile_line() */
#include <windows.h>
#endif

/* Reads a string of characters from fd while they match the string s.
   Returns the number of matching characters that were read. 
*/

static int 
facceptstr (FILE *fd, char *s)
{
  int c;
  unsigned i = 0;
  
  while (s[i] != 0 && (c = fgetc (fd)) != EOF)
    {
      if (c != s[i++])
        {
          ungetc (c, fd);
          return i-1;
        }
    }
  
  return i;
}

/* Accepts "\n" or "\r\n" or "\r". 
   Returns 1 if any of the three was read, 0 otherwise */

static int 
facceptnl (FILE *fd)
{
  int c, r = 0;

  c = fgetc (fd);
  if (c == '\r')
    {
      c = fgetc (fd);
      r = 1;
    }

  if (c == '\n')
    r = 1;
  else if (c != EOF)
    ungetc (c, fd);

  return r;
}

/* Reads a string from fd until the character "delim" or newline is seen, or 
   "len" characters have been written to "s" (including terminating null), 
   or EOF is reached. The "delim" and newline characters are left on the 
   stream.
   If s is NULL, characters are read from fd but not written anywhere.
   Returns the number of characters read.
*/

static int 
freadstrn (FILE *fd, char *s, char delim, unsigned int len)
{
  unsigned int i = 0;
  int c;
  
  while (i + 1 < len && (c = fgetc (fd)) != EOF)
    if (c == delim || IS_NEWLINE(c))
      {
        ungetc (c, fd);
        break;
      }
    else
      if (s != NULL)
        s[i++] = (char) c;
  
  if (i < len && s != NULL)
    s[i++] = 0;
  
  return i;
}

/* Reads an assignment from a save file. Return 1 if an assignment was
   successfully read, 0 if there are no more lines to read (at EOF) 
*/

int 
read_resumefile_line (int *method, mpz_t x, mpcandi_t *n, mpz_t sigma, mpz_t A, 
        mpz_t x0, double *b1, char *program, char *who, char *rtime, 
        char *comment, FILE *fd)
{
  int a, have_method, have_x, have_z, have_n, have_sigma, have_a, have_b1, 
      have_checksum, have_qx;
  unsigned int saved_checksum;
  char tag[16];
  mpz_t z;
  
  while (!feof (fd))
    {
      /* Ignore empty lines */
      if (facceptnl (fd))
        {
          continue;
        }
      
      /* Ignore lines beginning with '#'*/
      if (facceptstr (fd, "#"))
        {
          while (!facceptnl (fd) && !feof (fd))
            fgetc (fd);
          continue;
        }
      
      if (feof (fd))
        break;
      
      have_method = have_x = have_z = have_n = have_sigma = have_a = 
                    have_b1 = have_qx = have_checksum = 0;

      /* Set optional fields to zero */
      mpz_set_ui (sigma, 0);
      mpz_set_ui (A, 0);
      if (program != NULL)
        program[0] = 0;
      if (who != NULL)
        who[0] = 0;
      if (rtime != NULL)
        rtime[0] = 0;
      if (comment != NULL)
        comment[0] = 0;

      while (!facceptnl (fd) && !feof (fd))
        {
          freadstrn (fd, tag, '=', 16);
          
          if (!facceptstr (fd, "="))
            {
              printf ("Save file line has no equal sign after: %s\n", tag);
              goto error;
            }
          
          if (strcmp (tag, "METHOD") == 0)
            {
              if (facceptstr (fd, "ECM") == 3)
                {
                  *method = ECM_ECM;
                }
              else if (facceptstr (fd, "P"))
                {
                  a = facceptstr (fd, "-1");
                  if (a == 2)
                    {
                      *method = ECM_PM1;
                    }
                  else if (a == 0 && facceptstr (fd, "+1") == 2)
                    {
                      *method = ECM_PP1;
                    }
                  else 
                    goto error;
                }
              else
                goto error;

              have_method = 1;
            }
          else if (strcmp (tag, "X") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_x = 1;
            }
          else if (strcmp (tag, "Z") == 0)
            {
              mpz_init (z);
              mpz_inp_str (z, fd, 0);
              have_z = 1;
            }
          else if (strcmp (tag, "QX") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_qx = 1;
            }
          else if (strcmp (tag, "X0") == 0)
            {
              mpz_inp_str (x0, fd, 0);
            }
          else if (strcmp (tag, "CHECKSUM") == 0)
            {
              if (fscanf (fd, "%u", &saved_checksum) != 1)
                goto error;
              have_checksum = 1;
            }
          else if (strcmp (tag, "COMMENT") == 0)
            {
              freadstrn (fd, comment, ';', 255);
            }
          else if (strcmp (tag, "N") == 0)
            {
              /*mpz_inp_str (n, fd, 0);*/
	      /* we want to "maintain" any expressions, which were possibly stored in the file for N */
              have_n = read_number (n, fd, 0);
            }
          else if (strcmp (tag, "SIGMA") == 0)
            {
              mpz_inp_str (sigma, fd, 0);
              have_sigma = 1;
            }
          else if (strcmp (tag, "A") == 0)
            {
              mpz_inp_str (A, fd, 0);
              have_a = 1;
            }
          else if (strcmp (tag, "B1") == 0)
            {
              if (fscanf (fd, "%lf", b1) != 1)
                goto error;
              have_b1 = 1;
            }
          else if (strcmp (tag, "PROGRAM") == 0)
            {
              freadstrn (fd, program, ';', 255);
            }
          else if (strcmp (tag, "WHO") == 0)
            {
              freadstrn (fd, who, ';', 255);
            }
          else if (strcmp (tag, "TIME") == 0)
            {
              freadstrn (fd, rtime, ';', 255);
            }
          else /* Not a tag we know about */
            {
              printf ("Save file line has unknown tag: %s\n", tag);
              goto error;
            }
         
          /* Prime95 lines have no semicolon after SIGMA */
          if (!facceptstr (fd, ";") && ! (have_qx && have_n && have_sigma))
            {
              printf ("%s field not followed by semicolon\n", tag);
              goto error;
            }
          
          while (facceptstr (fd, " "));
        }
      
      /* Finished reading tags */
      
      /* Handle Prime95 v22 lines. These have no METHOD=ECM field and
         QX= instead of X= */
      
      if (have_qx)
        {
          if (have_n && have_sigma)
            {
              *method = ECM_ECM;
              /* *b1 = 1.0; */
              strcpy (program, "Prime95");
              mpz_mod (x, x, n->n);
              return 1;
            }
          goto error;
        }

#ifdef DEBUG
      if (*method != ECM_ECM && (have_sigma || have_a || have_z))
        {
          int count = have_sigma + have_a + have_z;
          printf ("Warning: Save file line has");
          if (have_sigma)
            {
              printf (" SIGMA");
              mpz_set_ui (sigma, 0);
              if (--count > 1)
                printf (",");
              else if (count > 0)
                printf (" and");
            }
          if (have_a)
            {
              printf (" A");
              mpz_set_ui (A, 0);
              if (--count > 0)
                printf (" and");
            }
          if (have_z)
            {
              printf (" Z");
              mpz_clear (Z);
              have_z = 0;
            }
          printf (" value for method other than ECM.\n");
        }
#endif
      
      if (!have_method || !have_x || !have_n || !have_b1 ||
          (method == ECM_ECM && !have_sigma && !have_a))
        {
          fprintf (stderr, "Save file line lacks fields\n");
          continue;
        }

      if (have_checksum)
        {
          mpz_t checksum;
          
          mpz_init (checksum);
          mpz_set_d (checksum, *b1);
          if (have_sigma)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
          if (have_a)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
          if (have_z)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (z, CHKSUMMOD));
          if (mpz_fdiv_ui (checksum, CHKSUMMOD) != saved_checksum)
            {
              fprintf (stderr, "Resume file line has bad checksum %u, expected %lu\n", 
                       saved_checksum, mpz_fdiv_ui (checksum, CHKSUMMOD));
              mpz_clear (checksum);
              continue;
            }
          mpz_clear (checksum);
        }

      mpz_mod (x, x, n->n);
      if (have_z)	/* Must normalize */
        {
          if (!mpz_invert (z, z, n->n)) /* Factor found? */
            {
              /* Oh great. What do we do with it now? */
              /* mpres_gcd (f, z, n); */
              printf ("Oops, factor found while reading from save file.\n");
            }
          mpz_mul (z, z, x);
          mpz_mod (x, z, n->n);
        }
        
      return 1;
      
error:
      /* In case of error, read rest of line and try next line */
      while (!facceptnl (fd) && !feof (fd))
        fgetc (fd);
    }
    
    /* We hit EOF without reading a proper save line */
    return 0;
}


/* Append a residue to the savefile with name given in fn.
   Returns 1 on success, 0 on error */
int  
write_resumefile_line (char *fn, int method, double B1, mpz_t sigma, mpz_t A, 
	mpz_t x, mpcandi_t *n, mpz_t x0, const char *comment)
{
  FILE *file;
  mpz_t checksum;
  time_t t;
  char text[256];
  char *uname, mname[32];
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  struct flock lock;
  int r, fd;
#endif

#ifdef DEBUG
  if (fn == NULL)
    {
      fprintf (stderr, "write_resumefile_line: fn == NULL\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  file = fopen (fn, "a");
  if (file == NULL)
    {
      fprintf (stderr, "Could not open file %s for writing\n", fn);
      return 0;
    }
  
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  /* Try to get a lock on the file so several processes can append to
     the same file safely */
  
  /* Supposedly some implementations of fcntl() can get confused over
     garbage in unused fields in a flock struct, so zero it */
  memset (&lock, 0, sizeof (struct flock));
  fd = fileno (file);
  lock.l_type = F_WRLCK;
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 1; 
  /* F_SETLKW: blocking exclusive lock request */
  r = fcntl (fd, F_SETLKW, &lock);
  if (r != 0)
    {
      fclose (file);
      return 0;
    }

  fseek (file, 0, SEEK_END);
#endif
  
  mpz_init (checksum);
  mpz_set_d (checksum, B1);
  fprintf (file, "METHOD=");
  if (method == ECM_PM1)
    fprintf (file, "P-1");
  else if (method == ECM_PP1)
    fprintf (file, "P+1");
  else 
    {
      fprintf (file, "ECM");
      if (mpz_sgn (sigma) != 0)
        {
          fprintf (file, "; SIGMA=");
          mpz_out_str (file, 10, sigma);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
        }
      else if (mpz_sgn (A) != 0)
        {
          fprintf (file, "; A=");
          mpz_out_str (file, 10, A);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
        }
    }
  
  fprintf (file, "; B1=%.0f; N=", B1);
  if (n->cpExpr)
    fprintf(file, "%s", n->cpExpr);
  else
    mpz_out_str (file, 10, n->n);
  fprintf (file, "; X=0x");
  mpz_out_str (file, 16, x);
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
#ifdef GPUECM
  fprintf (file, "; CHECKSUM=%lu; PROGRAM=GPU-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION_GPU);
#else
  fprintf (file, "; CHECKSUM=%lu; PROGRAM=GMP-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION);
#endif
  mpz_clear (checksum);
  
  if (mpz_sgn (x0) != 0)
    {
      fprintf (file, " X0=0x");
      mpz_out_str (file, 16, x0);
      fprintf (file, ";");
    }
  
  /* Try to get the users and his machines name */
  /* TODO: how to make portable? */
  uname = getenv ("LOGNAME");
  if (uname == NULL)
    uname = getenv ("USERNAME");
  if (uname == NULL)
    uname = "";
  
#if defined (_MSC_VER) || defined (__MINGW32__)
  /* dummy block, so that the vars needed here don't need to
    "spill" over to the rest of the function. */
  {
    DWORD size, i;
    TCHAR T[MAX_COMPUTERNAME_LENGTH+2];
    size=MAX_COMPUTERNAME_LENGTH+1;
    if (!GetComputerName(T, &size))
      strcpy(mname, "localPC");
    else
      {
        for (i = 0; i < sizeof(mname)-1; ++i)
          mname[i] = T[i];
        mname[sizeof(mname)-1] = 0;
      }
  }
#else
  if (gethostname (mname, 32) != 0)
    mname[0] = 0;
  mname[31] = 0; /* gethostname() may omit trailing 0 if hostname >31 chars */
#endif
  
  if (uname[0] != 0 || mname[0] != 0)
    {
      fprintf (file, " WHO=%.233s@%.32s;", uname, mname);
    }

  if (comment[0] != 0)
    fprintf (file, " COMMENT=%.255s;", comment);
  
  t = time (NULL);
  strncpy (text, ctime (&t), 255);
  text[255] = 0;
  text[strlen (text) - 1] = 0; /* Remove newline */
  fprintf (file, " TIME=%s;", text);
  fprintf (file, "\n");
  fflush (file);
#if defined(HAVE_FCNTL) && defined(HAVE_FILENO)
  lock.l_type = F_UNLCK;
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 1;  
  fcntl (fd, F_SETLKW, &lock); /* F_SETLKW: blocking lock request */
#endif
  fclose (file);

  return 1;
}
