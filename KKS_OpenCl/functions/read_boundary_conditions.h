#ifndef READ_BOUNDARY_CONDITIONS_H_
#define READ_BOUNDARY_CONDITIONS_H_

void read_boundary_conditions(char *argv[]) {
  FILE *fr;
  int i;
  char tempbuff[1000];
  
  char tmpstr1[1000];
  char tmpstr2[1000];
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
    fprintf(stderr, "Input file not found: %s\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  while(fgets(tempbuff,1000,fr)) {
    sscanf(tempbuff, "%999s = %999[^;];", tmpstr1, tmpstr2);
//     printf("%s\n",  tmpstr1);
//     printf("%s\n",  tmpstr2);
    if(tmpstr1[0] != '#') {
      if ((strcmp(tmpstr1, "BOUNDARY") == 0) && (NUMPHASES > 0)) {
        initialize_boundary_conditions(tmpstr2);
      }
      //else if ((strcmp(tmpstr1, "BOUNDARY_VALUE") == 0) && (NUMPHASES > 0)) {
        //initialize_boundary_points_values(tmpstr2);
      //}
    }
  }
  fclose(fr);
  
  char outfile[1004];
  
  if (snprintf(tmpstr2, sizeof(tmpstr2), "%s", argv[1])
      >= (int)sizeof(tmpstr2)) {
    fprintf(stderr, "Input filename is too long\n");
    exit(EXIT_FAILURE);
  }
  
  strcpy(tmpstr1,strtok(tmpstr2, "."));
  
  if (snprintf(outfile, sizeof(outfile), "%s.bd", tmpstr1)
      >= (int)sizeof(outfile)) {
    fprintf(stderr, "Derived boundary filename is too long\n");
    exit(EXIT_FAILURE);
  }
  
  fr = fopen(outfile, "w");
  if (fr == NULL) {
    fprintf(stderr, "Could not open derived boundary file: %s\n", outfile);
    exit(EXIT_FAILURE);
  }
  
  PRINT_BOUNDARY_CONDITIONS(fr);
  
  fclose(fr);
}
#endif
