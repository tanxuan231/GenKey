/*!
 ***********************************************************************
 *  \file
 *      mbuffer_mvc.c
 *
 *  \brief
 *      Frame buffer functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Athanasios Leontaris            <aleon@dolby.com>
 *      - Karsten Suehring
 *      - Alexis Tourapis                 <alexismt@ieee.org>
 *      - Yuwen He                        <yhe@dolby.com>
 ***********************************************************************
 */

#include <limits.h>

#include "global.h"
#include "header.h"
#include "image.h"
#include "mbuffer.h"
#include "mbuffer_common.h"
#include "mbuffer_mvc.h"
#include "memalloc.h"

#if (MVC_EXTENSION_ENABLE)

/*!
 ************************************************************************
 * \brief
 *    Reordering process for short-term reference pictures
 *
 ************************************************************************
 */
void reorder_short_term(Slice *currSlice, int cur_list, int num_ref_idx_lX_active_minus1, int picNumLX, int *refIdxLX, int currViewID)
{
  StorablePicture **RefPicListX = currSlice->listX[cur_list]; 
  int cIdx, nIdx;

  StorablePicture *picLX;

  picLX = get_short_term_pic(currSlice, currSlice->p_Dpb, picNumLX);

  for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
    if (RefPicListX[ cIdx ])
    {
      if( (RefPicListX[ cIdx ]->is_long_term ) ||  (RefPicListX[ cIdx ]->pic_num != picNumLX ) ||  
        ( currViewID != -1 && RefPicListX[ cIdx ]->layer_id  != currViewID )
        )
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];
    }
}


/*!
 ************************************************************************
 * \brief
 *    Reordering process for long-term reference pictures
 *
 ************************************************************************
 */
void reorder_long_term(Slice *currSlice, StorablePicture **RefPicListX, int num_ref_idx_lX_active_minus1, int LongTermPicNum, int *refIdxLX, int currViewID)
{
  int cIdx, nIdx;

  StorablePicture *picLX;

  picLX = get_long_term_pic(currSlice, currSlice->p_Dpb, LongTermPicNum);

  for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
    RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

  RefPicListX[ (*refIdxLX)++ ] = picLX;

  nIdx = *refIdxLX;

  for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
  {
    if (RefPicListX[ cIdx ])
    {
      if( (!RefPicListX[ cIdx ]->is_long_term ) ||  (RefPicListX[ cIdx ]->long_term_pic_num != LongTermPicNum ) ||  
          ( currViewID != -1 && RefPicListX[ cIdx ]->layer_id  != currViewID )
        )
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Returns inter-view prediction pic with given targetViewID
 *
 ************************************************************************
 */
static StorablePicture*  get_inter_view_pic(VideoParameters *p_Vid, Slice *currSlice, int targetViewID, int currPOC, int listidx)
{
  unsigned i;
  unsigned int listinterview_size;
  FrameStore **fs_listinterview;

  if (listidx == 0)
  {
    fs_listinterview = currSlice->fs_listinterview0;
    listinterview_size = currSlice->listinterviewidx0; 
  }
  else
  {
    fs_listinterview = currSlice->fs_listinterview1;
    listinterview_size = currSlice->listinterviewidx1; 
  }

  for(i=0; i<listinterview_size; i++)
  {
    if (fs_listinterview[i]->layer_id == GetVOIdx( p_Vid, targetViewID ))
    {
      if(p_Vid->structure==FRAME && fs_listinterview[i]->frame->poc == currPOC)
      {
        return fs_listinterview[i]->frame;
      }
      else if(p_Vid->structure==TOP_FIELD && fs_listinterview[i]->top_field->poc == currPOC)
      {
        return fs_listinterview[i]->top_field;
      }
      else if(p_Vid->structure==BOTTOM_FIELD && fs_listinterview[i]->bottom_field->poc == currPOC)
      {
        return fs_listinterview[i]->bottom_field;
      }
    }
  }

  return NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Generates a alternating field list from a given FrameStore inter-view list
 *
 ************************************************************************
 */
static void gen_pic_list_from_frame_interview_list(PictureStructure currStructure, FrameStore **fs_list, int list_idx, StorablePicture **list, char *list_size)
{
  int i;

  if (currStructure == TOP_FIELD)
  {
    for (i=0; i<list_idx; i++)
    {
      list[(int)(*list_size)] = fs_list[i]->top_field;
      (*list_size)++;
    }
  }
  if (currStructure == BOTTOM_FIELD)
  {
    for (i=0; i<list_idx; i++)
    {
      list[(int)(*list_size)] = fs_list[i]->bottom_field;
      (*list_size)++;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reordering process for inter-view reference pictures
 *
 ************************************************************************
 */
static void reorder_interview(VideoParameters *p_Vid, Slice *currSlice, StorablePicture **RefPicListX, int num_ref_idx_lX_active_minus1, int *refIdxLX, int targetViewID, int currPOC, int listidx)
{
  int cIdx, nIdx;
  StorablePicture *picLX;

  picLX = get_inter_view_pic(p_Vid, currSlice, targetViewID, currPOC, listidx);

  if (picLX)
  {
    for( cIdx = num_ref_idx_lX_active_minus1+1; cIdx > *refIdxLX; cIdx-- )
      RefPicListX[ cIdx ] = RefPicListX[ cIdx - 1];

    RefPicListX[ (*refIdxLX)++ ] = picLX;

    nIdx = *refIdxLX;

    for( cIdx = *refIdxLX; cIdx <= num_ref_idx_lX_active_minus1+1; cIdx++ )
    {
      if((GetViewIdx( p_Vid, RefPicListX[cIdx]->view_id ) != targetViewID) || (RefPicListX[cIdx]->poc != currPOC))
        RefPicListX[ nIdx++ ] = RefPicListX[ cIdx ];
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reordering process for MVC reference picture lists
 *
 ************************************************************************
 */
void reorder_ref_pic_list_mvc(Slice *currSlice, int cur_list, int **anchor_ref, int **non_anchor_ref, int view_id, int anchor_pic_flag, int currPOC, int listidx)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  int *modification_of_pic_nums_idc = currSlice->modification_of_pic_nums_idc[cur_list];
  int *abs_diff_pic_num_minus1 = currSlice->abs_diff_pic_num_minus1[cur_list];
  int *long_term_pic_idx = currSlice->long_term_pic_idx[cur_list];
  int num_ref_idx_lX_active_minus1 = currSlice->num_ref_idx_active[cur_list] - 1;
  int *abs_diff_view_idx_minus1 = currSlice->abs_diff_view_idx_minus1[cur_list];

  int i;

  int maxPicNum, currPicNum, picNumLXNoWrap, picNumLXPred, picNumLX;
  int picViewIdxLX, targetViewID;
  int refIdxLX = 0;
  int maxViewIdx =0;
  int curr_VOIdx = -1;
  int picViewIdxLXPred=-1;

  if (p_Vid->structure==FRAME)
  {
    maxPicNum  = p_Vid->max_frame_num;
    currPicNum = currSlice->frame_num;
  }
  else
  {
    maxPicNum  = 2 * p_Vid->max_frame_num;
    currPicNum = 2 * currSlice->frame_num + 1;
  }

  if(currSlice->svc_extension_flag==0)
  {
    curr_VOIdx = view_id;
    maxViewIdx = get_maxViewIdx(p_Vid, view_id, anchor_pic_flag, 0);
    picViewIdxLXPred=-1;
  }

  picNumLXPred = currPicNum;

  for (i=0; modification_of_pic_nums_idc[i]!=3; i++)
  {
    if (modification_of_pic_nums_idc[i] > 5)
      error ("Invalid modification_of_pic_nums_idc command", 500);

    if (modification_of_pic_nums_idc[i] < 2)
    {
      if (modification_of_pic_nums_idc[i] == 0)
      {
        if( picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 ) < 0 )
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 ) + maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred - ( abs_diff_pic_num_minus1[i] + 1 );
      }
      else // (modification_of_pic_nums_idc[i] == 1)
      {
        if( picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 )  >=  maxPicNum )
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 ) - maxPicNum;
        else
          picNumLXNoWrap = picNumLXPred + ( abs_diff_pic_num_minus1[i] + 1 );
      }
      picNumLXPred = picNumLXNoWrap;

      if( picNumLXNoWrap > currPicNum )
        picNumLX = picNumLXNoWrap - maxPicNum;
      else
        picNumLX = picNumLXNoWrap;

      reorder_short_term(currSlice, cur_list, num_ref_idx_lX_active_minus1, picNumLX, &refIdxLX, view_id);
    }
    else if (modification_of_pic_nums_idc[i] == 2) //(modification_of_pic_nums_idc[i] == 2)
    {
      reorder_long_term(currSlice, currSlice->listX[cur_list], num_ref_idx_lX_active_minus1, long_term_pic_idx[i], &refIdxLX, view_id);
    }
    else 
    {
      if(modification_of_pic_nums_idc[i] == 4) //(modification_of_pic_nums_idc[i] == 4)
      {
        picViewIdxLX = picViewIdxLXPred - (abs_diff_view_idx_minus1[i] + 1);
        if( picViewIdxLX <0)
          picViewIdxLX += maxViewIdx;
      }
      else //(modification_of_pic_nums_idc[i] == 5)
      {
        picViewIdxLX = picViewIdxLXPred + (abs_diff_view_idx_minus1[i] + 1);
        if( picViewIdxLX >= maxViewIdx)
          picViewIdxLX -= maxViewIdx;
      }
      picViewIdxLXPred = picViewIdxLX;

      if (anchor_pic_flag)
        targetViewID = anchor_ref[curr_VOIdx][picViewIdxLX];
      else
        targetViewID = non_anchor_ref[curr_VOIdx][picViewIdxLX];

      reorder_interview(p_Vid, currSlice, currSlice->listX[cur_list], num_ref_idx_lX_active_minus1, &refIdxLX, targetViewID, currPOC, listidx);
    }
  }
  // that's a definition
  currSlice->listXsize[cur_list] = (char) (num_ref_idx_lX_active_minus1 + 1);
}

void reorder_lists_mvc(Slice * currSlice, int currPOC)
{
  VideoParameters *p_Vid = currSlice->p_Vid;

  if ((currSlice->slice_type != I_SLICE)&&(currSlice->slice_type != SI_SLICE))
  {
    if (currSlice->ref_pic_list_reordering_flag[LIST_0])
    {
      reorder_ref_pic_list_mvc(currSlice, LIST_0,
        p_Vid->active_subset_sps->anchor_ref_l0,
        p_Vid->active_subset_sps->non_anchor_ref_l0,
        currSlice->view_id, currSlice->anchor_pic_flag, currPOC, 0);
    }
    if (p_Vid->no_reference_picture == currSlice->listX[0][currSlice->num_ref_idx_active[LIST_0]-1])
    {
      if (p_Vid->non_conforming_stream)
        printf("RefPicList0[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_active[LIST_0] - 1);
      else
        error("RefPicList0[ num_ref_idx_l0_active_minus1 ] in MVC layer is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    currSlice->listXsize[0] = (char)currSlice->num_ref_idx_active[LIST_0];
  }
  if (currSlice->slice_type == B_SLICE)
  {
    if (currSlice->ref_pic_list_reordering_flag[LIST_1])
    {
      reorder_ref_pic_list_mvc(currSlice, LIST_1,
        p_Vid->active_subset_sps->anchor_ref_l1,
        p_Vid->active_subset_sps->non_anchor_ref_l1,
        currSlice->view_id, currSlice->anchor_pic_flag, currPOC, 1);
    }
    if (p_Vid->no_reference_picture == currSlice->listX[1][currSlice->num_ref_idx_active[LIST_1]-1])
    {
      if (p_Vid->non_conforming_stream)
        printf("RefPicList1[ %d ] is equal to 'no reference picture'\n", currSlice->num_ref_idx_active[LIST_1] - 1);
      else
        error("RefPicList1[ num_ref_idx_l1_active_minus1 ] is equal to 'no reference picture', invalid bitstream",500);
    }
    // that's a definition
    currSlice->listXsize[1] = (char)currSlice->num_ref_idx_active[LIST_1];
  }

  free_ref_pic_list_reordering_buffer(currSlice);

  if ( currSlice->slice_type == P_SLICE )
  {    
#if PRINTREFLIST  
    unsigned int i;
    // print out for debug purpose
    if((p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) && currSlice->current_slice_nr==0)
    {
      if(currSlice->listXsize[0]>0)
      {
        printf("\n");
        printf(" ** (FinalViewID:%d) %s Ref Pic List 0 ****\n", currSlice->view_id, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
        for(i=0; i<(unsigned int)(currSlice->listXsize[0]); i++)  //ref list 0
        {
          printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[0][i]->poc, currSlice->listX[0][i]->pic_num, currSlice->listX[0][i]->view_id);
        }
      }
    }
#endif
  }
  else if ( currSlice->slice_type == B_SLICE )
  {
#if PRINTREFLIST
    unsigned int i;
    // print out for debug purpose
    if((p_Vid->profile_idc == MVC_HIGH || p_Vid->profile_idc == STEREO_HIGH) && currSlice->current_slice_nr==0)
    {
      if((currSlice->listXsize[0]>0) || (currSlice->listXsize[1]>0))
        printf("\n");
      if(currSlice->listXsize[0]>0)
      {
        printf(" ** (FinalViewID:%d) %s Ref Pic List 0 ****\n", currSlice->view_id, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
        for(i=0; i<(unsigned int)(currSlice->listXsize[0]); i++)  //ref list 0
        {
          printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[0][i]->poc, currSlice->listX[0][i]->pic_num, currSlice->listX[0][i]->view_id);
        }
      }
      if(currSlice->listXsize[1]>0)
      {
        printf(" ** (FinalViewID:%d) %s Ref Pic List 1 ****\n", currSlice->view_id, currSlice->structure==FRAME ? "FRM":(currSlice->structure==TOP_FIELD ? "TOP":"BOT"));
        for(i=0; i<(unsigned int)(currSlice->listXsize[1]); i++)  //ref list 1
        {
          printf("   %2d -> POC: %4d PicNum: %4d ViewID: %d\n", i, currSlice->listX[1][i]->poc, currSlice->listX[1][i]->pic_num, currSlice->listX[1][i]->view_id);
        }
      }
    }
#endif
  }
}

#endif
 
