
/*!
 *************************************************************************************
 * \file loopFilter.c
 *
 * \brief
 *    Filter to reduce blocking artifacts on a macroblock level.
 *    The filter strength is QP dependent.
 *
 * \author
 *    Contributors:
 *    - Peter List       Peter.List@t-systems.de:  Original code                                 (13-Aug-2001)
 *    - Jani Lainema     Jani.Lainema@nokia.com:   Some bug fixing, removal of recursiveness     (16-Aug-2001)
 *    - Peter List       Peter.List@t-systems.de:  inplace filtering and various simplifications (10-Jan-2002)
 *    - Anthony Joch     anthony@ubvideo.com:      Simplified switching between filters and
 *                                                 non-recursive default filter.                 (08-Jul-2002)
 *    - Cristina Gomila  cristina.gomila@thomson.net: Simplification of the chroma deblocking
 *                                                    from JVT-E089                              (21-Nov-2002)
 *    - Alexis Michael Tourapis atour@dolby.com:   Speed/Architecture improvements               (08-Feb-2007)
 *************************************************************************************
 */

#include "global.h"
#include "image.h"
#include "mb_access.h"

// likely already set - see testing via asserts
static void init_neighbors(VideoParameters *p_Vid)
{
  int i, j;
  int width = p_Vid->PicWidthInMbs;
  int height = p_Vid->PicHeightInMbs;
  int size = p_Vid->PicSizeInMbs;
  Macroblock *currMB = &p_Vid->mb_data[0];
  // do the top left corner
  currMB->mbup = NULL;
  currMB->mbleft = NULL;
  currMB++;
  // do top row
  for (i = 1; i < width; i++) 
  {
    currMB->mbup = NULL;
    currMB->mbleft = currMB - 1;
    currMB++;
  }

  // do left edge
  for (i = width; i < size; i += width) 
  {
    currMB->mbup = currMB - width;
    currMB->mbleft = NULL;   
    currMB += width;
  }
  // do all others
  for (j = width + 1; j < width * height + 1; j += width) 
  {
    currMB = &p_Vid->mb_data[j];
    for (i = 1; i < width; i++) 
    {
      currMB->mbup   = currMB - width;
      currMB->mbleft = currMB - 1;
      currMB++;
    }
  }
}


void  init_Deblock(VideoParameters *p_Vid, int mb_aff_frame_flag)
{
  if(p_Vid->yuv_format == YUV444 && p_Vid->separate_colour_plane_flag)
  {
    change_plane_JV(p_Vid, PLANE_Y, NULL);
    init_neighbors(p_Dec->p_Vid);
    change_plane_JV(p_Vid, PLANE_U, NULL);
    init_neighbors(p_Dec->p_Vid);
    change_plane_JV(p_Vid, PLANE_V, NULL);
    init_neighbors(p_Dec->p_Vid);
    change_plane_JV(p_Vid, PLANE_Y, NULL);
  }
  else 
    init_neighbors(p_Dec->p_Vid);	//ÔÝÊ±²»É¾

}

