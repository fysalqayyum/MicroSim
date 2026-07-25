#ifndef READ_BOUNDARY_CONDITIONS_H_
#define READ_BOUNDARY_CONDITIONS_H_

void read_boundary_conditions(char *argv[]) {
  FILE *fr;
  int i;
  char tempbuff[1000];
  
  char tmpstr1[100];
  char tmpstr2[100];
  char **tmp;
  
  bool decision;
  
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  
  long k, j;
  long index;
  long length;
  long phase;
  
  fr = fopen(argv[1], "rt");
  if(fr == NULL) {
    printf("\nERROR!!! Stop. \n- [%d] Infile `%s` not found!!\n\n",taskid, argv[1]);
    exit(1);
  }
  if (taskid == MASTER){
    printf("\n- [%d] Reading boundary conditions from `%s`",taskid, argv[1]);
  }
  
  
  while(fgets(tempbuff,1000,fr)) {
    if (sscanf(tempbuff, "%99s = %99[^;];", tmpstr1, tmpstr2) != 2) {
      continue;
    }
//     printf("%s\n",  tmpstr1);
//     printf("%s\n",  tmpstr2);
    if(tmpstr1[0] != '#') {
      if ((strcmp(tmpstr1, "BOUNDARY") == 0) && (NUMPHASES > 0)) {
        initialize_boundary_conditions(tmpstr2);
      }
      else if ((strcmp(tmpstr1, "BOUNDARY_VALUE") == 0) && (NUMPHASES > 0)) {
        initialize_boundary_points_values(tmpstr2);
      }
    }
  }
  fclose(fr);
}
#endif
