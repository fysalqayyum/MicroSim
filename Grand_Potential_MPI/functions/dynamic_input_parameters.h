#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>

#include "utility_functions_new.h"

typedef struct microsimDynamicInputFileParameters_t
{
  char    file_path[1000];
  time_t  last_init_time;

} msDyInp_t;

msDyInp_t ms_dyInp;

void msDyInp_Initialize(char* argv[])
{ 
  strcpy(ms_dyInp.file_path, argv[1]);

  if(taskid == MASTER){
    printf("\n- [%d] Dynamic Input File set to be : %s", taskid, ms_dyInp.file_path) ;
  }

  struct stat file_stat;
  if (stat(ms_dyInp.file_path, &file_stat) != 0) {
      perror("stat");
  }
  ms_dyInp.last_init_time = file_stat.st_mtime;
}

bool msDyInp_IsModified()
{
    struct stat file_stat;
    if (stat(ms_dyInp.file_path, &file_stat) != 0) {
        perror("stat");
        return false;
    }

    time_t last_modified = file_stat.st_mtime;

    if (ms_dyInp.last_init_time != last_modified) {
      printf("\n- [%d] Input file %s was recently modified.", taskid, ms_dyInp.file_path);
      ms_dyInp.last_init_time = last_modified ;
      return true ;
    } else {
      return false ;
    }
}


void msDyInp_ReInitializeInpParameters()
{
  FILE * File_ptr = ms_FILE_Open(ms_dyInp.file_path, "rt");

  printf("\n- [%d] re initializing variables",taskid);
  
  char tempbuff[10000];
  char tmpstr1[10000];
  char tmpstr2[10000];
  
  while(fgets(tempbuff, 10000, File_ptr))
  {
    if (sscanf(tempbuff, "%9999s = %9999[^;];", tmpstr1, tmpstr2) != 2) {
      continue;
    }

    if (tmpstr1[0] != '#')
    {
      if (strcmp(tmpstr1,"DELTA_t")==0) {
        deltat = atof(tmpstr2);
      }
      else if (strcmp(tmpstr1,"NTIMESTEPS")==0) {
        ntimesteps = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"SAVET")==0) {
        saveT = atol(tmpstr2);
      }
      else if (strcmp(tmpstr1,"TRACK_PROGRESS")==0) {
        time_output = atol(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Rotation_matrix") == 0) && (FUNCTION_ANISOTROPY))
      {
        populate_rotation_matrix(Rotation_matrix, Inv_Rotation_matrix, tmpstr2);
      }
      else if ((strcmp(tmpstr1, "Tempgrady") == 0) && (TEMPGRADY))
      {
        char **tmp = (char**)malloc(sizeof(char*)*5);
        char *str1, *saveptr1, *token;
        for (i = 0; i < 5; ++i) {
          tmp[i] = (char*)malloc(sizeof(char)*10);
        }

        for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL) {
          token = strtok_r(str1, "{,}", &saveptr1);
          if (token == NULL)
              break;
          strcpy(tmp[i],token);
        }
        // temperature_gradientY.base_temp          = atof(tmp[0]);
        // temperature_gradientY.DeltaT             = atof(tmp[1]);
        // temperature_gradientY.Distance           = atof(tmp[2]);
        // temperature_gradientY.gradient_OFFSET    = atof(tmp[3]);
        temperature_gradientY.velocity           = atof(tmp[4]);
        for (i = 0; i < 5; ++i) {
          free(tmp[i]);
        }
        free(tmp);
      }
      else if ((strcmp(tmpstr1, "dab") == 0)) {
        populate_matrix(dab, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "fab") == 0)) {
        populate_matrix(fab, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "Gamma_abc") == 0)) {
        populate_matrix3M(Gamma_abc, tmpstr2, NUMPHASES);
      }
      else if ((strcmp(tmpstr1, "DIFFUSIVITY") == 0)) {
        populate_diffusivity_matrix(Diffusivity, tmpstr2, NUMCOMPONENTS);
      }
      else if (LBM && (strcmp(tmpstr1, "LBM_SAVE_FREQ") == 0)){
        lbmSaveFreq = atoi(tmpstr2);
      }
      else if (LBM && (strcmp(tmpstr1, "dt") == 0)) {
        dt = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "eidt_mode") == 0) && EIDT_SWITCH)
      {
        EIDT.mode = atol(tmpstr2);  
      }
      // else if ((strcmp(tmpstr1, "comp_ff") == 0) && EIDT_SWITCH)
      // {
      //   populate_double_array(EIDT.comp_far_field, tmpstr2, NUMCOMPONENTS-1);        
      // }
      else if ((strcmp(tmpstr1, "comp_ff_rate") == 0) && EIDT_SWITCH)
      {
        populate_double_array(EIDT.comp_rate, tmpstr2, NUMCOMPONENTS-1);
      }
      else if ((strcmp(tmpstr1, "rad_coeff") == 0) && EIDT_SWITCH)
      {
        EIDT.radius_coeff = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "step_coeff") == 0) && EIDT_SWITCH)
      {
        EIDT.step_coeff = atof(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "step_zero") == 0) && EIDT_SWITCH)
      {
        EIDT.step0 = atol(tmpstr2);
      }
      else {
        continue;
      }
    }
  }

  fclose(File_ptr);
}
