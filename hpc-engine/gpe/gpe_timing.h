/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   WARSAW UNIVERSITY OF TECHNOLOGY                                       *
 *   FACULTY OF PHYSICS                                                    *
 *   NUCLEAR THEORY GROUP                                                  *
 *   See also AUTHORS file                                                 *
 *                                                                         *
 *   This file is a part of GPE for GPU project.                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/ 

#include <sys/time.h>
static double t_gettimeofday ;
static clock_t t_clock , t_clock1 ;
static struct timeval s ;

void b_t( void ) 
{ /* hack together a clock w/ microsecond resolution */
  gettimeofday( &s , NULL ) ;
  t_clock = clock() ;
  t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec ;
}

double e_t( int type = 0 ) 
{
  switch ( type ) 
    {
    case 0 :
      t_clock1 = clock() ;
      gettimeofday( &s , NULL ) ;
      t_clock = t_clock1 - t_clock ;
      t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec - t_gettimeofday ;
      return t_gettimeofday ;
    case 1 :
      return t_gettimeofday ;
    case 2 :
      return t_clock / ( double ) CLOCKS_PER_SEC ;
    }
  return t_gettimeofday ;
}
