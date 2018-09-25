#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include "readpar.h"

ino_t LoadParFile (char *file);

struct RP_parfile_par {
  ino_t inode;
  char  *name;
  char  *value;
};

#define NPARMAX 100000
struct RP_parfile_par RP_parinfo[NPARMAX];
int RP_nread_par=0;

struct stat RP_sbuff_read[100];
int RP_nfiles_loaded=0;

#define LINELEN 10000

int get_parameter_file_par (char *file, char *parname, char *parvalue)
{
  int f;
  int lastarg, lastval;
  f = get_parameter_file_par_work (file,parname,parvalue);
  if (RP_q_expand_vars())
    if (defined(parvalue))
      RP_expand_variables(parvalue);
  return f;
}

int get_parameter_file_par_work (char *file, char *parname, char *parvalue) {
  ino_t inode;
  char name[200];
  char par[200];
  char *p;
  int qcase;
  int i;
  
  if ( ! (inode=LoadParFile (file)) ) {
    strcpy (parvalue,RP_UNDEFINED);
    return 0;
  }
  
  if (parname[0]=='^') {
    qcase = 1;
    strcpy(name,parname+1);
  }
  else {
    qcase = 0;
    strcpy (name,parname);
    lcase (name);
  }

  for (i=0;i<RP_nread_par;i++) {
    if ( RP_parinfo[i].inode == inode ) {
      if (!qcase) {
	strcpy(par,RP_parinfo[i].name);
	lcase(par); p=par;
      }
      else
	p = RP_parinfo[i].name;
      
      if (strcmp(p,name)==0) {
	strcpy(parvalue,RP_parinfo[i].value);
	return 1;
      }
    }
  }
  strcpy(parvalue,RP_UNDEFINED);
  return 0;
}

/*************************************************************************/
ino_t LoadParFile (char *file) {
  struct stat sbuff;
  char buff[200];
  char name[100],value[10000],comment[100],value0[10000];
  int valuecont;
  int i,k,l;
  FILE *parfile;
  ino_t inode;
  char line[LINELEN]; int type;
  char groupname[80], work[80], *p;
  int nold, nnew;

  /*  fprintf (stderr,"Reading par file %s\n",file); */

  /* first, stat the file */
  if (stat (file,&sbuff)) { /* if cannot be stat'ed, return */
    strcpy(buff,"Cannot read parameter file ");
    strcat(buff,file);
    perror(buff);
    return 0;
  }
  
  /* see if this inode was already read in */
  for (i=0;i<RP_nfiles_loaded;i++) {
    if ( sbuff.st_ino == RP_sbuff_read[i].st_ino ) { 
      /* already read in; was it modified since then? */
      if ( sbuff.st_ctime <= RP_sbuff_read[i].st_ctime ) {
	/* file was not modified */
	return sbuff.st_ino;
      }
      else {
	/* FREE UP OLD INFO */
	RP_Rem_PF_Pars (sbuff.st_ino);
	/* and immediate exit from this loop */
	i=RP_nfiles_loaded+1;
      }
    }
  }
  
  /* Open parameter file */
  parfile = fopen(file,"r");
  if ( parfile == NULL ) {
    strcpy(buff,"Cannot read parameter file ");
    strcat(buff,file);
    perror(buff);
    return 0;
  }
  
  /* Copy inode info */
  RP_nfiles_loaded++; k=RP_nfiles_loaded-1; RP_sbuff_read[k]=sbuff;
  inode = sbuff.st_ino;
  
  /* Read file */
  RP_IsInGroup(0);
  RP_IsPAR_Continued(0);
  
  while (fgets(line,LINELEN,parfile) != NULL) {
    type = RP_RecordType (line);

    switch (type) {

    case RP_BEGGROUP:  /* Beggining of new data group */
      RP_IsInGroup(1);
      RP_Parse_Group_String (line,groupname,comment);
      break;

    case RP_INCLSTRING: /* Include another par file */
      strcpy(buff,line+9); /* Include format is @include{filename} */
      i=0;
      while (buff[i]!='}'&&buff[i]!='\0') i++;
      buff[i]='\0';
      if (RP_q_expand_vars()) 
        RP_expand_variables(buff);
      nold = RP_nread_par;
      LoadParFile(buff); /* Recursive call of LoadParFile */
      nnew = RP_nread_par;
      /* now we need to fix inode in freshly read in parameters */
      for (i=nold;i<nnew;i++) {
	RP_parinfo[i].inode = inode;
      }
      break;

    case RP_ENDGROUP:  /* End of new data group */
      RP_IsInGroup(0);
      break;

    case RP_PARSTRING: /* Parameter string */
      RP_Parse_Par_String (line,name,value,comment,&valuecont);
      if ( valuecont ) {
	RP_IsPAR_Continued(1);
	strcpy (value0,value);
      }
      else
	RP_IsPAR_Continued(0);
      
      if (RP_PAR_Continued()) break;
      /* otherwise, store par and value */
      if (RP_InGroup()) { /* set par name to group.parname */
	strcpy (work,groupname); strcat (work,".");
	strcat (work,name); strcpy(name,work);
      }
      RP_Push_New_PF_Par(inode,name,value);
      break;

    case RP_PARCONT:  /* Continuation of parameter value */
      RP_Parse_Cont_String (line,value,comment,&valuecont);
      strcat (value0,"\n");
      strcat (value0,value);
      if ( valuecont )
	RP_IsPAR_Continued(1);
      else
	RP_IsPAR_Continued(0);
      
      if (RP_PAR_Continued()) break;
      /* otherwise, store par and value */
      if (RP_InGroup()) { /* set par name to group.parname */
	strcpy (work,groupname); strcat (work,".");
	strcat (work,name); strcpy(name,work);
      }
      RP_Push_New_PF_Par(inode,name,value0);
      break;

    case RP_IRAFSTRING:
      RP_Parse_IRAF_String (line,name,value,comment);
      RP_IsPAR_Continued(0);
      if (RP_InGroup()) { /* set par name to group.parname */
	strcpy (work,groupname); strcat (work,".");
	strcat (work,name); strcpy(name,work);
      }
      RP_Push_New_PF_Par(inode,name,value);
      break;
    }
  }
  fclose(parfile);
  

  /*  
      fprintf (stderr,"Par file %s loaded: %d entries total\n",file,RP_nread_par); 
  */


  return inode;
}  

int  RP_Push_New_PF_Par (ino_t inode, char *name, char *value)
{
  int l;
  RP_nread_par++; l=RP_nread_par-1;
  RP_parinfo[l].inode = inode;
  if ((RP_parinfo[l].name = 
       (char *) malloc(sizeof(char)*(strlen(name)+1)))==NULL){
    perror("Can't allocate par name space");
    exit(1);
  }
  strcpy(RP_parinfo[l].name,name);
  if ((RP_parinfo[l].value = 
       (char *) malloc(sizeof(char)*(strlen(value)+1)))==NULL){
    perror("Can't allocate par value space");
    exit(1);
  }
  strcpy(RP_parinfo[l].value,value);
}


int RP_Rem_PF_Pars (ino_t inode)
{
  int n,i;
  
  n=0;
  for (i=0;i<RP_nread_par;i++) {
    if ( RP_parinfo[i].inode != inode ) {
      if ( n < i ) {
	RP_parinfo[n]=RP_parinfo[i];
	n++;
      }
    }
    else { /* entry to remove */
      free (RP_parinfo[i].name);
      free (RP_parinfo[i].value);
    }
  }
  RP_nread_par=n;

  n=0;
  for (i=0;i<RP_nfiles_loaded;i++) {
    if ( inode != RP_sbuff_read[i].st_ino ) {
      RP_sbuff_read[n]=RP_sbuff_read[i];
      n++;
    }
  }
  RP_nfiles_loaded=n;
  /*  fprintf (stderr,"Removed inode %d; %d files loaded, %d parameters known\n",
      inode,RP_nfiles_loaded,RP_nread_par); */
}

int RP_Line_By_Line_Parse (char *line){	    
  /* Call RP_IsInGroup(0); RP_IsPAR_Continued(0); before iterating this prog */
  
  int type;
  static int valuecont;
  char name[80],value[1024],comment[1024];

  type = RP_RecordType (line);

  switch (type) { /* Change status setting if necessary */

  case RP_BEGGROUP:  /* Beggining of new data group */
      RP_IsInGroup(1);
      break;

  case RP_ENDGROUP:  /* End of new data group */
    RP_IsInGroup(0);
    break;

  case RP_PARSTRING: /* Parameter string */
    RP_Parse_Par_String (line,name,value,comment,&valuecont);
    if ( valuecont )
      RP_IsPAR_Continued(1);
    else
      RP_IsPAR_Continued(0);
      
    break;

  case RP_PARCONT:  /* Continuation of parameter value */
    RP_Parse_Cont_String (line,value,comment,&valuecont);
    if ( valuecont )
      RP_IsPAR_Continued(1);
    else
      RP_IsPAR_Continued(0);
    
    break;

  case RP_IRAFSTRING:
    RP_IsPAR_Continued(0);
    break;
  }
  return type;
}
