
/*!
 ***********************************************************************
 *  \file
 *      mbuffer.c
 *
 *  \brief
 *      Frame buffer functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Suehring
 *      - Alexis Tourapis                 <alexismt@ieee.org>
 *      - Jill Boyce                      <jill.boyce@thomson.net>
 *      - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *      - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *      - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *      - Yuwen He                        <yhe@dolby.com>
 ***********************************************************************
 */

#include <limits.h>

#include "global.h"
#include "erc_api.h"
#include "image.h"
#include "mbuffer.h"
#include "mbuffer_common.h"
#include "memalloc.h"
#include "fast_memory.h"
#include "input.h"


#define MAX_LIST_SIZE 33

void alloc_pic_motion(PicMotionParamsOld *motion, int size_y, int size_x)
{
  motion->mb_field = calloc (size_y * size_x, sizeof(byte));
  if (motion->mb_field == NULL)
    no_mem_exit("alloc_storable_picture: motion->mb_field");
}

/*!
 ************************************************************************
 * \brief
 *    Allocate memory for a stored picture.
 *
 * \param p_Vid
 *    VideoParameters
 * \param structure
 *    picture structure
 * \param size_x
 *    horizontal luma size
 * \param size_y
 *    vertical luma size
 * \param size_x_cr
 *    horizontal chroma size
 * \param size_y_cr
 *    vertical chroma size
 *
 * \return
 *    the allocated StorablePicture structure
 ************************************************************************
 */
StorablePicture* alloc_storable_picture(VideoParameters *p_Vid, PictureStructure structure, int size_x, int size_y, int size_x_cr, int size_y_cr, int is_output)
{
  seq_parameter_set_rbsp_t *active_sps = p_Vid->active_sps;  

  StorablePicture *s;
  int   nplane;

  //printf ("Allocating (%s) picture (x=%d, y=%d, x_cr=%d, y_cr=%d)\n", (type == FRAME)?"FRAME":(type == TOP_FIELD)?"TOP_FIELD":"BOTTOM_FIELD", size_x, size_y, size_x_cr, size_y_cr);

  s = calloc (1, sizeof(StorablePicture));
  if (NULL==s)
    no_mem_exit("alloc_storable_picture: s");

  if (structure!=FRAME)
  {
    size_y    /= 2;
    size_y_cr /= 2;
  }

  s->PicSizeInMbs = (size_x*size_y)/256;
  s->imgUV = NULL;

  get_mem2Dpel_pad (&(s->imgY), size_y, size_x, p_Vid->iLumaPadY, p_Vid->iLumaPadX);
  s->iLumaStride = size_x+2*p_Vid->iLumaPadX;
  s->iLumaExpandedHeight = size_y+2*p_Vid->iLumaPadY;

  if (active_sps->chroma_format_idc != YUV400)
  {
    get_mem3Dpel_pad(&(s->imgUV), 2, size_y_cr, size_x_cr, p_Vid->iChromaPadY, p_Vid->iChromaPadX);
  }

  s->iChromaStride =size_x_cr + 2*p_Vid->iChromaPadX;
  s->iChromaExpandedHeight = size_y_cr + 2*p_Vid->iChromaPadY;
  s->iLumaPadY   = p_Vid->iLumaPadY;
  s->iLumaPadX   = p_Vid->iLumaPadX;
  s->iChromaPadY = p_Vid->iChromaPadY;
  s->iChromaPadX = p_Vid->iChromaPadX;

  s->separate_colour_plane_flag = p_Vid->separate_colour_plane_flag;

  get_mem2Dmp     ( &s->mv_info, (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
  alloc_pic_motion( &s->motion , (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));

  if( (p_Vid->separate_colour_plane_flag != 0) )
  {
    for( nplane=0; nplane<MAX_PLANE; nplane++ )
    {
      get_mem2Dmp      (&s->JVmv_info[nplane], (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
      alloc_pic_motion(&s->JVmotion[nplane] , (size_y >> BLOCK_SHIFT), (size_x >> BLOCK_SHIFT));
    }
  }

  s->pic_num   = 0;
  s->frame_num = 0;
  s->long_term_frame_idx = 0;
  s->long_term_pic_num   = 0;
  s->used_for_reference  = 0;
  s->is_long_term        = 0;
  s->non_existing        = 0;
  s->is_output           = 0;
  s->max_slice_id        = 0;
#if (MVC_EXTENSION_ENABLE)
  s->view_id = -1;
#endif

  s->structure=structure;

  s->size_x = size_x;
  s->size_y = size_y;
  s->size_x_cr = size_x_cr;
  s->size_y_cr = size_y_cr;
  s->size_x_m1 = size_x - 1;
  s->size_y_m1 = size_y - 1;
  s->size_x_cr_m1 = size_x_cr - 1;
  s->size_y_cr_m1 = size_y_cr - 1;

  s->top_field    = p_Vid->no_reference_picture;
  s->bottom_field = p_Vid->no_reference_picture;
  s->frame        = p_Vid->no_reference_picture;

  s->dec_ref_pic_marking_buffer = NULL;

  s->coded_frame  = 0;
  s->mb_aff_frame_flag  = 0;

  s->top_poc = s->bottom_poc = s->poc = 0;
  s->seiHasTone_mapping = 0;

  if(!p_Vid->active_sps->frame_mbs_only_flag && structure != FRAME)
  {
    int i, j;
    for(j = 0; j < MAX_NUM_SLICES; j++)
    {
      for (i = 0; i < 2; i++)
      {
        s->listX[j][i] = calloc(MAX_LIST_SIZE, sizeof (StorablePicture*)); // +1 for reordering
        if (NULL==s->listX[j][i])
        no_mem_exit("alloc_storable_picture: s->listX[i]");
      }
    }
  }

  return s;
}

/*!
 ************************************************************************
 * \brief
 *    Free frame store memory.
 *
 * \param p_Vid
 *    VideoParameters
 * \param f
 *    FrameStore to be freed
 *
 ************************************************************************
 */
void free_frame_store(FrameStore* f)
{
  if (f)
  {
    if (f->frame)
    {
      free_storable_picture(f->frame);
      f->frame=NULL;
    }
    if (f->top_field)
    {
      free_storable_picture(f->top_field);
      f->top_field=NULL;
    }
    if (f->bottom_field)
    {
      free_storable_picture(f->bottom_field);
      f->bottom_field=NULL;
    }
    free(f);
  }
}

void free_pic_motion(PicMotionParamsOld *motion)
{
  if (motion->mb_field)
  {
    free(motion->mb_field);
    motion->mb_field = NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Free picture memory.
 *
 * \param p
 *    Picture to be freed
 *
 ************************************************************************
 */
void free_storable_picture(StorablePicture* p)
{
  int nplane;
  if (p)
  {
    if (p->mv_info)
    {
      free_mem2Dmp(p->mv_info);
      p->mv_info = NULL;
    }
    free_pic_motion(&p->motion);

    if( (p->separate_colour_plane_flag != 0) )
    {
      for( nplane=0; nplane<MAX_PLANE; nplane++ )
      {
        if (p->JVmv_info[nplane])
        {
          free_mem2Dmp(p->JVmv_info[nplane]);
          p->JVmv_info[nplane] = NULL;
        }
        free_pic_motion(&p->JVmotion[nplane]);
      }
    }

    if (p->imgY)
    {
      free_mem2Dpel_pad(p->imgY, p->iLumaPadY, p->iLumaPadX);
      p->imgY = NULL;
    }

    if (p->imgUV)
    {
      free_mem3Dpel_pad(p->imgUV, 2, p->iChromaPadY, p->iChromaPadX);
      p->imgUV=NULL;
    }


    if (p->seiHasTone_mapping)
      free(p->tone_mapping_lut);

    {
      int i, j;
      for(j = 0; j < MAX_NUM_SLICES; j++)
      {
        for(i=0; i<2; i++)
        {
          if(p->listX[j][i])
          {
            free(p->listX[j][i]);
            p->listX[j][i] = NULL;
          }
        }
      }
    }
    free(p);
    p = NULL;
  }
}
/*!
 ************************************************************************
 * \brief
 *    Initialize reference lists depending on current slice type
 *
 ************************************************************************
 */
void init_lists_i_slice(Slice *currSlice)
{

#if (MVC_EXTENSION_ENABLE)
  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;
#endif

  currSlice->listXsize[0] = 0;
  currSlice->listXsize[1] = 0;
}
/*!
 ************************************************************************
 * \brief
 *    Initialize reference lists for a P Slice
 *
 ************************************************************************
 */
void init_lists_p_slice(Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  DecodedPictureBuffer *p_Dpb = currSlice->p_Dpb;

  unsigned int i;

  int list0idx = 0;
  int listltidx = 0;

  FrameStore **fs_list0;
  FrameStore **fs_listlt;

#if (MVC_EXTENSION_ENABLE)
  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;
#endif

  if (currSlice->structure == FRAME)
  {
    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_used==3)
      {
        if ((p_Dpb->fs_ref[i]->frame->used_for_reference) && (!p_Dpb->fs_ref[i]->frame->is_long_term))
        {
          currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
        }
      }
    }
    // order list 0 by PicNum
    qsort((void *)currSlice->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_pic_num_desc);
    currSlice->listXsize[0] = (char) list0idx;
    //printf("listX[0] (PicNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", currSlice->listX[0][i]->pic_num);} printf("\n");

    // long term handling
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ltref[i]->is_used==3)
      {
        if (p_Dpb->fs_ltref[i]->frame->is_long_term)
        {
          currSlice->listX[0][list0idx++] = p_Dpb->fs_ltref[i]->frame;
        }
      }
    }
    qsort((void *)&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
    currSlice->listXsize[0] = (char) list0idx;
  }
  else
  {
    fs_list0 = calloc(p_Dpb->size, sizeof (FrameStore*));
    if (NULL==fs_list0)
      no_mem_exit("init_lists: fs_list0");
    fs_listlt = calloc(p_Dpb->size, sizeof (FrameStore*));
    if (NULL==fs_listlt)
      no_mem_exit("init_lists: fs_listlt");

    for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
    {
      if (p_Dpb->fs_ref[i]->is_reference)
      {
        fs_list0[list0idx++] = p_Dpb->fs_ref[i];
      }
    }

    qsort((void *)fs_list0, list0idx, sizeof(FrameStore*), compare_fs_by_frame_num_desc);

    //printf("fs_list0 (FrameNum): "); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->frame_num_wrap);} printf("\n");

    currSlice->listXsize[0] = 0;
    gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);

    //printf("listX[0] (PicNum): "); for (i=0; i < currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->pic_num);} printf("\n");

    // long term handling
    for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
    {
      fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
    }

    qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

    gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);

    free(fs_list0);
    free(fs_listlt);
  }
  currSlice->listXsize[1] = 0;


  // set max size
  currSlice->listXsize[0] = (char) imin (currSlice->listXsize[0], currSlice->num_ref_idx_active[LIST_0]);
  currSlice->listXsize[1] = (char) imin (currSlice->listXsize[1], currSlice->num_ref_idx_active[LIST_1]);

  // set the unused list entries to NULL
  for (i=currSlice->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[0][i] = p_Vid->no_reference_picture;
  }
  for (i=currSlice->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[1][i] = p_Vid->no_reference_picture;
  }

#if PRINTREFLIST
#if (MVC_EXTENSION_ENABLE)
  // print out for debug purpose
  if((p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) && currSlice->current_slice_nr==0)
  {
    if(currSlice->listXsize[0]>0)
    {
      printf("\n");
      printf(" ** (CurViewID:%d %d) %s Ref Pic List 0 ****\n", currSlice->view_id, currSlice->ThisPOC, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
      for(i=0; i<(unsigned int)(currSlice->listXsize[0]); i++)  //ref list 0
      {
        printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[0][i]->poc, currSlice->listX[0][i]->pic_num, currSlice->listX[0][i]->view_id);
      }
    }
  }
#endif
#endif
}


/*!
 ************************************************************************
 * \brief
 *    Initialize reference lists for a B Slice
 *
 ************************************************************************
 */
void init_lists_b_slice(Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  DecodedPictureBuffer *p_Dpb = currSlice->p_Dpb;

  unsigned int i;
  int j;

  int list0idx = 0;
  int list0idx_1 = 0;
  int listltidx = 0;

  FrameStore **fs_list0;
  FrameStore **fs_list1;
  FrameStore **fs_listlt;

#if (MVC_EXTENSION_ENABLE)
  currSlice->listinterviewidx0 = 0;
  currSlice->listinterviewidx1 = 0;
#endif

  {
    // B-Slice
    if (currSlice->structure == FRAME)
    {
      for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (currSlice->framepoc >= p_Dpb->fs_ref[i]->frame->poc) //!KS use >= for error concealment
            {
              currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)currSlice->listX[0], list0idx, sizeof(StorablePicture*), compare_pic_by_poc_desc);

      //get the backward reference picture (POC>current POC) in list0;
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used==3)
        {
          if ((p_Dpb->fs_ref[i]->frame->used_for_reference)&&(!p_Dpb->fs_ref[i]->frame->is_long_term))
          {
            if (currSlice->framepoc < p_Dpb->fs_ref[i]->frame->poc)
            {
              currSlice->listX[0][list0idx++] = p_Dpb->fs_ref[i]->frame;
            }
          }
        }
      }
      qsort((void *)&currSlice->listX[0][list0idx_1], list0idx-list0idx_1, sizeof(StorablePicture*), compare_pic_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        currSlice->listX[1][list0idx-list0idx_1+j]=currSlice->listX[0][j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        currSlice->listX[1][j-list0idx_1]=currSlice->listX[0][j];
      }

      currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;

      //printf("listX[0] (PicNum): "); for (i=0; i<currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->pic_num);} printf("\n");
      //printf("listX[1] (PicNum): "); for (i=0; i<currSlice->listXsize[1]; i++){printf ("%d  ", currSlice->listX[1][i]->pic_num);} printf("\n");
      //printf("currSlice->listX[0] currPoc=%d (Poc): ", p_Vid->framepoc); for (i=0; i<currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->poc);} printf("\n");
      //printf("currSlice->listX[1] currPoc=%d (Poc): ", p_Vid->framepoc); for (i=0; i<currSlice->listXsize[1]; i++){printf ("%d  ", currSlice->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ltref[i]->is_used==3)
        {
          if (p_Dpb->fs_ltref[i]->frame->is_long_term)
          {
            currSlice->listX[0][list0idx]   = p_Dpb->fs_ltref[i]->frame;
            currSlice->listX[1][list0idx++] = p_Dpb->fs_ltref[i]->frame;
          }
        }
      }
      qsort((void *)&currSlice->listX[0][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      qsort((void *)&currSlice->listX[1][(short) currSlice->listXsize[0]], list0idx - currSlice->listXsize[0], sizeof(StorablePicture*), compare_pic_by_lt_pic_num_asc);
      currSlice->listXsize[0] = currSlice->listXsize[1] = (char) list0idx;
    }
    else
    {
      fs_list0 = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_list0)
        no_mem_exit("init_lists: fs_list0");
      fs_list1 = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_list1)
        no_mem_exit("init_lists: fs_list1");
      fs_listlt = calloc(p_Dpb->size, sizeof (FrameStore*));
      if (NULL==fs_listlt)
        no_mem_exit("init_lists: fs_listlt");

      currSlice->listXsize[0] = 0;
      currSlice->listXsize[1] = 1;

      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (currSlice->ThisPOC >= p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)fs_list0, list0idx, sizeof(FrameStore*), compare_fs_by_poc_desc);
      list0idx_1 = list0idx;
      for (i=0; i<p_Dpb->ref_frames_in_buffer; i++)
      {
        if (p_Dpb->fs_ref[i]->is_used)
        {
          if (currSlice->ThisPOC < p_Dpb->fs_ref[i]->poc)
          {
            fs_list0[list0idx++] = p_Dpb->fs_ref[i];
          }
        }
      }
      qsort((void *)&fs_list0[list0idx_1], list0idx-list0idx_1, sizeof(FrameStore*), compare_fs_by_poc_asc);

      for (j=0; j<list0idx_1; j++)
      {
        fs_list1[list0idx-list0idx_1+j]=fs_list0[j];
      }
      for (j=list0idx_1; j<list0idx; j++)
      {
        fs_list1[j-list0idx_1]=fs_list0[j];
      }

      //printf("fs_list0 currPoc=%d (Poc): ", currSlice->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list0[i]->poc);} printf("\n");
      //printf("fs_list1 currPoc=%d (Poc): ", currSlice->ThisPOC); for (i=0; i<list0idx; i++){printf ("%d  ", fs_list1[i]->poc);} printf("\n");

      currSlice->listXsize[0] = 0;
      currSlice->listXsize[1] = 0;
      gen_pic_list_from_frame_list(currSlice->structure, fs_list0, list0idx, currSlice->listX[0], &currSlice->listXsize[0], 0);
      gen_pic_list_from_frame_list(currSlice->structure, fs_list1, list0idx, currSlice->listX[1], &currSlice->listXsize[1], 0);

      //printf("currSlice->listX[0] currPoc=%d (Poc): ", p_Vid->framepoc); for (i=0; i<currSlice->listXsize[0]; i++){printf ("%d  ", currSlice->listX[0][i]->poc);} printf("\n");
      //printf("currSlice->listX[1] currPoc=%d (Poc): ", p_Vid->framepoc); for (i=0; i<currSlice->listXsize[1]; i++){printf ("%d  ", currSlice->listX[1][i]->poc);} printf("\n");

      // long term handling
      for (i=0; i<p_Dpb->ltref_frames_in_buffer; i++)
      {
        fs_listlt[listltidx++]=p_Dpb->fs_ltref[i];
      }

      qsort((void *)fs_listlt, listltidx, sizeof(FrameStore*), compare_fs_by_lt_pic_idx_asc);

      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, currSlice->listX[0], &currSlice->listXsize[0], 1);
      gen_pic_list_from_frame_list(currSlice->structure, fs_listlt, listltidx, currSlice->listX[1], &currSlice->listXsize[1], 1);

      free(fs_list0);
      free(fs_list1);
      free(fs_listlt);
    }
  }

  if ((currSlice->listXsize[0] == currSlice->listXsize[1]) && (currSlice->listXsize[0] > 1))
  {
    // check if lists are identical, if yes swap first two elements of currSlice->listX[1]
    int diff=0;
    for (j = 0; j< currSlice->listXsize[0]; j++)
    {
      if (currSlice->listX[0][j] != currSlice->listX[1][j])
      {
        diff = 1;
        break;
      }
    }
    if (!diff)
    {
      StorablePicture *tmp_s = currSlice->listX[1][0];
      currSlice->listX[1][0]=currSlice->listX[1][1];
      currSlice->listX[1][1]=tmp_s;
    }
  }

  // set max size
  currSlice->listXsize[0] = (char) imin (currSlice->listXsize[0], currSlice->num_ref_idx_active[LIST_0]);
  currSlice->listXsize[1] = (char) imin (currSlice->listXsize[1], currSlice->num_ref_idx_active[LIST_1]);

  // set the unused list entries to NULL
  for (i=currSlice->listXsize[0]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[0][i] = p_Vid->no_reference_picture;
  }
  for (i=currSlice->listXsize[1]; i< (MAX_LIST_SIZE) ; i++)
  {
    currSlice->listX[1][i] = p_Vid->no_reference_picture;
  }

#if PRINTREFLIST
#if (MVC_EXTENSION_ENABLE)
  // print out for debug purpose
  if((p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) && currSlice->current_slice_nr==0)
  {
    if((currSlice->listXsize[0]>0) || (currSlice->listXsize[1]>0))
      printf("\n");
    if(currSlice->listXsize[0]>0)
    {
      printf(" ** (CurViewID:%d %d) %s Ref Pic List 0 ****\n", currSlice->view_id, currSlice->ThisPOC, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
      for(i=0; i<(unsigned int)(currSlice->listXsize[0]); i++)  //ref list 0
      {
        printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[0][i]->poc, currSlice->listX[0][i]->pic_num, currSlice->listX[0][i]->view_id);
      }
    }
    if(currSlice->listXsize[1]>0)
    {
      printf(" ** (CurViewID:%d %d) %s Ref Pic List 1 ****\n", currSlice->view_id, currSlice->ThisPOC, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
      for(i=0; i<(unsigned int)(currSlice->listXsize[1]); i++)  //ref list 1
      {
        printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[1][i]->poc, currSlice->listX[1][i]->pic_num, currSlice->listX[1][i]->view_id);
      }
    }
  }
#endif
#endif
}
 /*!
 ************************************************************************
 * \brief
 *    Returns short term pic with given picNum
 *
 ************************************************************************
 */
StorablePicture*  get_short_term_pic(Slice *currSlice, DecodedPictureBuffer *p_Dpb, int picNum)
{
  unsigned i;

  for (i = 0; i < p_Dpb->ref_frames_in_buffer; i++)
  {
    if (currSlice->structure == FRAME)
    {
      if (p_Dpb->fs_ref[i]->is_reference == 3)
        if ((!p_Dpb->fs_ref[i]->frame->is_long_term)&&(p_Dpb->fs_ref[i]->frame->pic_num == picNum))
          return p_Dpb->fs_ref[i]->frame;
    }
    else
    {
      if (p_Dpb->fs_ref[i]->is_reference & 1)
        if ((!p_Dpb->fs_ref[i]->top_field->is_long_term)&&(p_Dpb->fs_ref[i]->top_field->pic_num == picNum))
          return p_Dpb->fs_ref[i]->top_field;
      if (p_Dpb->fs_ref[i]->is_reference & 2)
        if ((!p_Dpb->fs_ref[i]->bottom_field->is_long_term)&&(p_Dpb->fs_ref[i]->bottom_field->pic_num == picNum))
          return p_Dpb->fs_ref[i]->bottom_field;
    }
  }

  return currSlice->p_Vid->no_reference_picture;
}
/*!
 ************************************************************************
 * \brief
 *    Allocate memory for buffering of reference picture reordering commands
 ************************************************************************
 */
void alloc_ref_pic_list_reordering_buffer(Slice *currSlice)
{
  if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE)
  {
    int size = currSlice->num_ref_idx_active[LIST_0] + 1;
    if ((currSlice->modification_of_pic_nums_idc[LIST_0] = calloc(size ,sizeof(int)))==NULL) 
       no_mem_exit("alloc_ref_pic_list_reordering_buffer: modification_of_pic_nums_idc_l0");
    if ((currSlice->abs_diff_pic_num_minus1[LIST_0] = calloc(size,sizeof(int)))==NULL) 
       no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l0");
    if ((currSlice->long_term_pic_idx[LIST_0] = calloc(size,sizeof(int)))==NULL) 
       no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l0");
#if (MVC_EXTENSION_ENABLE)
    if ((currSlice->abs_diff_view_idx_minus1[LIST_0] = calloc(size,sizeof(int)))==NULL) 
       no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_view_idx_minus1_l0");
#endif
  }
  else
  {
    currSlice->modification_of_pic_nums_idc[LIST_0] = NULL;
    currSlice->abs_diff_pic_num_minus1[LIST_0] = NULL;
    currSlice->long_term_pic_idx[LIST_0] = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->abs_diff_view_idx_minus1[LIST_0] = NULL;
#endif
  }

  if (currSlice->slice_type == B_SLICE)
  {
    int size = currSlice->num_ref_idx_active[LIST_1] + 1;
    if ((currSlice->modification_of_pic_nums_idc[LIST_1] = calloc(size,sizeof(int)))==NULL) 
      no_mem_exit("alloc_ref_pic_list_reordering_buffer: modification_of_pic_nums_idc_l1");
    if ((currSlice->abs_diff_pic_num_minus1[LIST_1] = calloc(size,sizeof(int)))==NULL) 
      no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l1");
    if ((currSlice->long_term_pic_idx[LIST_1] = calloc(size,sizeof(int)))==NULL) 
      no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l1");
#if (MVC_EXTENSION_ENABLE)
    if ((currSlice->abs_diff_view_idx_minus1[LIST_1] = calloc(size,sizeof(int)))==NULL) 
      no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_view_idx_minus1_l1");
#endif
  }
  else
  {
    currSlice->modification_of_pic_nums_idc[LIST_1] = NULL;
    currSlice->abs_diff_pic_num_minus1[LIST_1] = NULL;
    currSlice->long_term_pic_idx[LIST_1] = NULL;
#if (MVC_EXTENSION_ENABLE)
    currSlice->abs_diff_view_idx_minus1[LIST_1] = NULL;
#endif
  }
}
/*!
 ************************************************************************
 * \brief
 *    Free memory for buffering of reference picture reordering commands
 ************************************************************************
 */
void free_ref_pic_list_reordering_buffer(Slice *currSlice)
{
  if (currSlice->modification_of_pic_nums_idc[LIST_0])
    free(currSlice->modification_of_pic_nums_idc[LIST_0]);
  if (currSlice->abs_diff_pic_num_minus1[LIST_0])
    free(currSlice->abs_diff_pic_num_minus1[LIST_0]);
  if (currSlice->long_term_pic_idx[LIST_0])
    free(currSlice->long_term_pic_idx[LIST_0]);

  currSlice->modification_of_pic_nums_idc[LIST_0] = NULL;
  currSlice->abs_diff_pic_num_minus1[LIST_0] = NULL;
  currSlice->long_term_pic_idx[LIST_0] = NULL;

  if (currSlice->modification_of_pic_nums_idc[LIST_1])
    free(currSlice->modification_of_pic_nums_idc[LIST_1]);
  if (currSlice->abs_diff_pic_num_minus1[LIST_1])
    free(currSlice->abs_diff_pic_num_minus1[LIST_1]);
  if (currSlice->long_term_pic_idx[LIST_1])
    free(currSlice->long_term_pic_idx[LIST_1]);

  currSlice->modification_of_pic_nums_idc[LIST_1] = NULL;
  currSlice->abs_diff_pic_num_minus1[LIST_1] = NULL;
  currSlice->long_term_pic_idx[LIST_1] = NULL;

#if (MVC_EXTENSION_ENABLE)
  if (currSlice->abs_diff_view_idx_minus1[LIST_0])
    free(currSlice->abs_diff_view_idx_minus1[LIST_0]);
  currSlice->abs_diff_view_idx_minus1[LIST_0] = NULL;
  if (currSlice->abs_diff_view_idx_minus1[LIST_1])
    free(currSlice->abs_diff_view_idx_minus1[LIST_1]);
  currSlice->abs_diff_view_idx_minus1[LIST_1] = NULL;
#endif
}
static int is_view_id_in_ref_view_list(int view_id, int *ref_view_id, int num_ref_views)
{
   int i;
   for(i=0; i<num_ref_views; i++)
   {
     if(view_id == ref_view_id[i])
       break;
   }

   return (num_ref_views && (i<num_ref_views));
}
void append_interview_list(DecodedPictureBuffer *p_Dpb, 
                           PictureStructure currPicStructure, //0: frame; 1:top field; 2: bottom field;
                           int list_idx, 
                           FrameStore **list,
                           int *listXsize, 
                           int currPOC, 
                           int curr_view_id, 
                           int anchor_pic_flag)
{
  VideoParameters *p_Vid = p_Dpb->p_Vid;
  int iVOIdx = curr_view_id;
  int pic_avail;
  int poc = 0;
  int fld_idx;
  int num_ref_views, *ref_view_id;
  FrameStore *fs = p_Dpb->fs_ilref[0];


  if(iVOIdx <0)
    printf("Error: iVOIdx: %d is not less than 0\n", iVOIdx);

  if(anchor_pic_flag)
  {
    num_ref_views = list_idx? p_Vid->active_subset_sps->num_anchor_refs_l1[iVOIdx] : p_Vid->active_subset_sps->num_anchor_refs_l0[iVOIdx];
    ref_view_id   = list_idx? p_Vid->active_subset_sps->anchor_ref_l1[iVOIdx]:p_Vid->active_subset_sps->anchor_ref_l0[iVOIdx];
  }
  else
  {
    num_ref_views = list_idx? p_Vid->active_subset_sps->num_non_anchor_refs_l1[iVOIdx] : p_Vid->active_subset_sps->num_non_anchor_refs_l0[iVOIdx];
    ref_view_id = list_idx? p_Vid->active_subset_sps->non_anchor_ref_l1[iVOIdx]:p_Vid->active_subset_sps->non_anchor_ref_l0[iVOIdx];
  }

  //  if(num_ref_views <= 0)
  //    printf("Error: iNumOfRefViews: %d is not larger than 0\n", num_ref_views);

  if(currPicStructure == BOTTOM_FIELD)
    fld_idx = 1;
  else
    fld_idx = 0;

    if(currPicStructure==FRAME)
    {
      pic_avail = (fs->is_used == 3);
      if (pic_avail)
        poc = fs->frame->poc;
    }
    else if(currPicStructure==TOP_FIELD)
    {
      pic_avail = fs->is_used & 1;
      if (pic_avail)
        poc = fs->top_field->poc;
    }
    else if(currPicStructure==BOTTOM_FIELD)
    {
      pic_avail = fs->is_used & 2;
      if (pic_avail)
        poc = fs->bottom_field->poc;
    }
    else
      pic_avail =0;

    if(pic_avail && fs->inter_view_flag[fld_idx])
    {
      if(poc == currPOC)
      {
        if(is_view_id_in_ref_view_list(fs->view_id, ref_view_id, num_ref_views))
        {
          //add one inter-view reference;
          list[*listXsize] = fs; 
          //next;
          (*listXsize)++;
        }
      }
    }
}
